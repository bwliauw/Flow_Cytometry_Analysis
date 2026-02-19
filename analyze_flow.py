import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import os
import argparse
import base64
from io import BytesIO
import numpy as np

REQUIRED_METRIC_KEYS = ("percent_parent", "mfi_ratio", "mfi_af488")
REQUIRED_MAPPING_COLUMNS = ("Sample Name", "Updated Sample Name", "Sample Type")


def clean_and_merge(data_dir):
    # 1. Locate files
    files = os.listdir(data_dir)

    def pick_csvs_from_dir(dir_path):
        local_raw = None
        local_mapping = None
        for name in os.listdir(dir_path):
            if name.endswith(".csv"):
                if "FlowJo table" in name and local_raw is None:
                    local_raw = os.path.join(dir_path, name)
                elif "plate_mapping" in name and local_mapping is None:
                    local_mapping = os.path.join(dir_path, name)
        return local_raw, local_mapping

    raw_csv, mapping_csv = pick_csvs_from_dir(data_dir)
    if not raw_csv or not mapping_csv:
        for child in files:
            child_path = os.path.join(data_dir, child)
            if os.path.isdir(child_path):
                nested_raw, nested_mapping = pick_csvs_from_dir(child_path)
                if nested_raw and nested_mapping:
                    raw_csv, mapping_csv = nested_raw, nested_mapping
                    break
    
    if not raw_csv or not mapping_csv:
        raise FileNotFoundError(f"Could not find required CSV files in {data_dir}")
    
    print(f"Loading raw data from: {raw_csv}")
    print(f"Loading mapping from: {mapping_csv}")
    
    # 2. Load and clean raw data
    # The first column is sample names, remaining are outputs
    raw_df = pd.read_csv(raw_csv)
    
    # Remove "Mean" and "SD" rows (last two rows usually, but let's be safe)
    # The user said: "The last two rows contain artifact values of 'mean' and 'SD'"
    # We can filter by checking if the first column contains "Mean" or "SD"
    sample_col = raw_df.columns[0]
    normalized_sample_names = raw_df[sample_col].astype(str).str.strip().str.lower()
    raw_df = raw_df[~normalized_sample_names.isin(["mean", "sd"])]
    # Also remove any trailing empty rows if they exist
    raw_df = raw_df.dropna(subset=[sample_col])
    
    # 3. Load mapping data
    mapping_df = pd.read_csv(mapping_csv)
    mapping_col_map = resolve_mapping_columns(mapping_df.columns)

    # Rename resolved mapping columns to canonical names so downstream code is stable.
    mapping_df = mapping_df.rename(
        columns={
            mapping_col_map["Sample Name"]: "Sample Name",
            mapping_col_map["Updated Sample Name"]: "Updated Sample Name",
            mapping_col_map["Sample Type"]: "Sample Type",
        }
    )
    
    # 4. Merge
    # We'll merge on the first column of raw_df and 'Sample Name' in mapping_df
    merged_df = pd.merge(
        raw_df, 
        mapping_df, 
        left_on=sample_col, 
        right_on='Sample Name', 
        how='left'
    )

    # Fail fast if any raw sample was not mapped; otherwise groupby would silently drop NaN keys.
    unmatched_mask = merged_df["Updated Sample Name"].isna() | merged_df["Sample Type"].isna()
    if unmatched_mask.any():
        unmatched_samples = (
            merged_df.loc[unmatched_mask, sample_col]
            .astype(str)
            .drop_duplicates()
            .tolist()
        )
        preview = ", ".join(unmatched_samples[:10])
        suffix = "..." if len(unmatched_samples) > 10 else ""
        raise ValueError(
            "Found raw samples missing mapping metadata (Updated Sample Name and/or Sample Type). "
            f"Unmatched samples ({len(unmatched_samples)}): {preview}{suffix}"
        )
    
    # Save intermediate CSV
    output_path = "processed_flow_data.csv"
    merged_df.to_csv(output_path, index=False)
    print(f"Saved merged data to {output_path}")
    
    return merged_df

def identify_columns(df):
    cols = df.columns
    target_cols = {}
    
    # 1) Singlets/AF647(+)/AF488(+) %Parent
    # Pattern: Freq. of Parent (%)
    parent_match = [c for c in cols if 'Freq. of Parent (%)' in str(c) and 'AF488(+)' in str(c)]
    if parent_match:
        target_cols['percent_parent'] = parent_match[-1] # Usually the more specific one
    
    # 2) Singlets/AF647(+)/AF488(+) MFI_(AF488/AF647) single-cell ratio
    # Pattern: Ratio_AF488_AF647
    ratio_match = [c for c in cols if 'Ratio_AF488_AF647' in str(c)]
    if ratio_match:
        target_cols['mfi_ratio'] = ratio_match[-1]
        
    # 3) Singlets/AF647(+)/AF488(+) MFI_AF488
    # Pattern: AF488-A
    af488_match = [c for c in cols if 'AF488-A' in str(c) and 'AF488(+)' in str(c)]
    if af488_match:
        target_cols['mfi_af488'] = af488_match[-1]
        
    return target_cols


def _normalize_column_name(name):
    return "".join(ch for ch in str(name).strip().lower() if ch.isalnum())


def resolve_mapping_columns(mapping_columns):
    normalized = {_normalize_column_name(col): col for col in mapping_columns}
    required_norm = {_normalize_column_name(col): col for col in REQUIRED_MAPPING_COLUMNS}

    resolved = {}
    missing = []
    for norm_required, canonical in required_norm.items():
        matched_col = normalized.get(norm_required)
        if matched_col is None:
            missing.append(canonical)
        else:
            resolved[canonical] = matched_col

    if missing:
        available_preview = ", ".join([str(col) for col in mapping_columns])
        raise ValueError(
            "Mapping CSV is missing required columns. "
            f"Missing: {', '.join(missing)}. "
            f"Available columns: {available_preview}"
        )

    return resolved


def validate_target_columns(target_cols, available_columns):
    missing = [key for key in REQUIRED_METRIC_KEYS if key not in target_cols]
    if not missing:
        return

    missing_desc = ", ".join(missing)
    available_preview = ", ".join([str(col) for col in available_columns])
    raise ValueError(
        "Unable to identify all required metric columns. "
        f"Missing keys: {missing_desc}. "
        "Please verify the input CSV has expected FlowJo output columns. "
        f"Available columns: {available_preview}"
    )

def get_sem(x):
    return x.std() / np.sqrt(len(x)) if len(x) > 1 else 0

def generate_plots(df, target_cols, export_png=False):
    validate_target_columns(target_cols, df.columns)

    # Color mapping
    color_map = {
        'Negative Control': '#D3D3D3', # Light Grey
        'Positive Control': '#FFB6C1', # Light Red
        'Experimental Sample': '#ADD8E6' # Light Blue
    }
    
    # Group by Sample Name and Type for plotting averages
    # We want to keep the order from the mapping file if possible
    # We use 'Updated Sample Name' and 'Sample Type' from the mapping file
    plot_data = df.groupby(['Updated Sample Name', 'Sample Type'], sort=False).agg({
        target_cols['percent_parent']: ['mean', get_sem],
        target_cols['mfi_ratio']: ['mean', get_sem],
        target_cols['mfi_af488']: ['mean', get_sem]
    }).reset_index()
    
    # Flatten columns
    plot_data.columns = ['Sample Name', 'Sample Type', 
                         'percent_parent_mean', 'percent_parent_sem',
                         'mfi_ratio_mean', 'mfi_ratio_sem',
                         'mfi_af488_mean', 'mfi_af488_sem']
    
    metrics = [
        ('percent_parent', 'Singlets/AF647(+)/AF488(+) %Parent'),
        ('mfi_ratio', 'Singlets/AF647(+)/AF488(+) MFI_(AF488/AF647) single-cell ratio'),
        ('mfi_af488', 'Singlets/AF647(+)/AF488(+) MFI_AF488')
    ]
    
    plot_base64 = []
    
    for metric_id, title in metrics:
        fig, ax = plt.subplots(figsize=(14, 7))
        x_positions = np.arange(len(plot_data))
        means = plot_data[f"{metric_id}_mean"].to_numpy()
        sems = plot_data[f"{metric_id}_sem"].to_numpy()
        bar_colors = [color_map.get(sample_type, "#B0B0B0") for sample_type in plot_data["Sample Type"]]

        # Draw bars and error bars together to keep a strict 1:1 row-to-bar mapping.
        ax.bar(
            x_positions,
            means,
            yerr=sems,
            capsize=3,
            color=bar_colors,
            edgecolor="black",
            linewidth=0.5,
        )

        ax.set_title(title)
        ax.set_xticks(x_positions)
        ax.set_xticklabels(plot_data["Sample Name"], rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Value (Mean ± SEM)")
        ax.set_xlabel("Sample Name")

        legend_order = ["Negative Control", "Positive Control", "Experimental Sample"]
        present_types = plot_data["Sample Type"].dropna().astype(str).unique().tolist()
        legend_handles = [
            Patch(facecolor=color_map[sample_type], edgecolor="black", label=sample_type)
            for sample_type in legend_order
            if sample_type in present_types
        ]
        if legend_handles:
            ax.legend(handles=legend_handles, title="Sample Type", bbox_to_anchor=(1.05, 1), loc="upper left")
        fig.tight_layout()
        
        # Save to base64
        buf = BytesIO()
        fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
        buf.seek(0)
        img_str = base64.b64encode(buf.read()).decode('utf-8')
        plot_base64.append(img_str)
        
        if export_png:
            filename = f"{metric_id}_plot.png"
            fig.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Exported {filename}")
            
        plt.close(fig)
        
    return plot_base64, plot_data

def generate_report(plot_data, plot_base64, target_cols):
    with open("experiment_summary.md", "w") as f:
        f.write("# Flow Cytometry Analysis Summary\n\n")
        
        f.write("## Data Table\n\n")
        # Format table for markdown
        table_df = plot_data.copy()
        
        # Map back to original column names for the table headers if desired, 
        # but let's use cleaner names for the table.
        display_cols = {
            'Sample Name': 'Sample Name',
            'Sample Type': 'Sample Type',
            'percent_parent': 'Singlets/AF647(+)/AF488(+) %Parent',
            'mfi_ratio': 'MFI Ratio (AF488/AF647)',
            'mfi_af488': 'MFI AF488'
        }
        
        for metric in ['percent_parent', 'mfi_ratio', 'mfi_af488']:
            table_df[display_cols[metric]] = table_df.apply(lambda r: f"{r[metric+'_mean']:.2f} ± {r[metric+'_sem']:.2f}", axis=1)
        
        final_table = table_df[['Sample Name', 'Sample Type'] + [display_cols[m] for m in ['percent_parent', 'mfi_ratio', 'mfi_af488']]]
        f.write(final_table.to_markdown(index=False))
        f.write("\n\n")
        
        f.write("## Figures\n\n")
        titles = [
            "Singlets/AF647(+)/AF488(+) %Parent",
            "MFI Ratio (AF488/AF647)",
            "MFI AF488"
        ]
        for title, img in zip(titles, plot_base64):
            f.write(f"### {title}\n")
            f.write(f"![{title}](data:image/png;base64,{img})\n\n")
            
    print("Generated experiment_summary.md")

def main():
    parser = argparse.ArgumentParser(description="Analyze Flow Cytometry Data")
    parser.add_argument("data_dir", help="Path to directory containing raw CSV and plate mapping CSV")
    parser.add_argument("--export-png", action="store_true", help="Export plots as PNG files")
    args = parser.parse_args()
    
    try:
        merged_df = clean_and_merge(args.data_dir)
        target_cols = identify_columns(merged_df)
        validate_target_columns(target_cols, merged_df.columns)
        
        print("Identified target columns:")
        for k, v in target_cols.items():
            print(f"  {k}: {v}")
            
        plot_base64, plot_data = generate_plots(merged_df, target_cols, args.export_png)
        generate_report(plot_data, plot_base64, target_cols)
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
