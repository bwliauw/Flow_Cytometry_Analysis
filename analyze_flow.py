"""Flow cytometry analysis pipeline.

This script takes a directory containing two CSV inputs:
1) Raw FlowJo-exported results.
2) Plate mapping metadata (sample name mapping, sample type, replicate).

It produces three output artifacts in a project-local output folder:
- `processed_flow_data.csv` (cleaned + merged per-row data)
- `experiment_summary.md` (key findings + summary table + figure links)
- Four PNG figures (bar charts with mean ± SEM)

Core responsibilities:
- Locate and validate expected input files.
- Clean known raw-data artifacts ("Mean"/"SD" rows).
- Resolve mapping column names robustly even with formatting variations.
- Fail fast when mapping metadata is incomplete.
- Compute requested thresholds and summary metrics.
- Render consistently ordered figures with styling and threshold overlays.
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import os
import argparse
import numpy as np

REQUIRED_METRIC_KEYS = ("percent_parent", "mfi_ratio", "mfi_af488", "mirfp_expression")
REQUIRED_MAPPING_COLUMNS = ("Sample Name", "Updated Sample Name", "Sample Type", "Replicate")


def get_output_dir(data_dir):
    """Return/create the project-local output directory for a given dataset path.

    The output folder name is derived from the input folder name:
    `<input_folder>_analyzed_data`
    and is always created inside the same directory as this script.
    """
    input_dir = os.path.abspath(os.path.normpath(data_dir))
    input_name = os.path.basename(input_dir)
    project_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(project_dir, f"{input_name}_analyzed_data")
    os.makedirs(output_dir, exist_ok=True)
    return output_dir


def clean_and_merge(data_dir, output_dir):
    """Load raw + mapping CSVs, clean artifacts, merge metadata, and export CSV.

    Steps:
    - Discover raw and mapping CSVs (top-level first, then one nested level).
    - Remove known artifact rows ("Mean", "SD") from raw data.
    - Normalize/validate required mapping columns.
    - Left-merge mapping metadata onto each raw sample row.
    - Fail explicitly if any sample lacks required metadata.
    - Reformat columns for the exported intermediate CSV.
    """
    # 1. Locate files
    files = os.listdir(data_dir)

    def pick_csvs_from_dir(dir_path):
        """Pick the first matching raw and mapping CSV from a directory.

        Matching is intentionally flexible for raw file names to handle
        naming differences across exports (FlowJo table / ratio reanalysis).
        """
        local_raw = None
        local_mapping = None
        for name in os.listdir(dir_path):
            if name.endswith(".csv"):
                lower_name = name.lower()
                if "plate_mapping" in lower_name and local_mapping is None:
                    local_mapping = os.path.join(dir_path, name)
                elif local_raw is None and (
                    "flowjo table" in lower_name
                    or "ratio_reanalysis" in lower_name
                    or "reanaly" in lower_name
                ):
                    local_raw = os.path.join(dir_path, name)
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
            mapping_col_map["Replicate"]: "Replicate",
        }
    )
    
    # 4. Merge
    # We'll merge on the first column of raw_df and 'Sample Name' in mapping_df
    # Left join preserves every raw row so analysis never silently drops
    # a sample that appeared in the instrument export.
    merged_df = pd.merge(
        raw_df, 
        mapping_df, 
        left_on=sample_col, 
        right_on='Sample Name', 
        how='left'
    )

    # Fail fast if any raw sample was not mapped; otherwise groupby would silently drop NaN keys.
    unmatched_mask = (
        merged_df["Updated Sample Name"].isna()
        | merged_df["Sample Type"].isna()
        | merged_df["Replicate"].isna()
    )
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
            "Found raw samples missing mapping metadata (Updated Sample Name, Sample Type, and/or Replicate). "
            f"Unmatched samples ({len(unmatched_samples)}): {preview}{suffix}"
        )

    # Keep raw sample IDs in the first column and remove duplicated merge key from mapping table.
    if sample_col != "Sample Name":
        merged_df = merged_df.drop(columns=["Sample Name"], errors="ignore")
        merged_df = merged_df.rename(columns={sample_col: "Sample Name"})

    # Rename column for user-facing clarity.
    merged_df = merged_df.rename(columns={"Updated Sample Name": "True Sample Name"})

    # Remove unused export columns from FlowJo and preserve meaningful data columns.
    drop_unnamed_cols = [col for col in merged_df.columns if str(col).startswith("Unnamed:")]
    if drop_unnamed_cols:
        merged_df = merged_df.drop(columns=drop_unnamed_cols, errors="ignore")

    # Ensure metadata columns are front-loaded in the output CSV.
    # This keeps the exported table human-readable and consistent between runs.
    leading_cols = ["Sample Name", "True Sample Name", "Sample Type", "Replicate"]
    remaining_cols = [col for col in merged_df.columns if col not in leading_cols]
    merged_df = merged_df[leading_cols + remaining_cols]
    
    # Save intermediate CSV
    output_path = os.path.join(output_dir, "processed_flow_data.csv")
    merged_df.to_csv(output_path, index=False)
    print(f"Saved merged data to {output_path}")
    
    return merged_df

def identify_columns(df):
    """Identify FlowJo metric columns using pattern-based matching.

    Returns:
        dict with canonical keys:
          - percent_parent
          - mfi_ratio
          - mfi_af488
          - mirfp_expression

    The script uses canonical keys downstream so analysis logic is independent
    from exact raw CSV column ordering.
    """
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

    # 4) Singlets/AF647(+) Geometric Mean (R1-A :: miRFP-A)
    mirfp_match = [
        c for c in cols
        if "a-FLAG_AF647(+)" in str(c)
        and "Geometric Mean (R1-A :: miRFP-A)" in str(c)
        and "/a-His_AF488(+)" not in str(c)
    ]
    if mirfp_match:
        target_cols["mirfp_expression"] = mirfp_match[-1]
        
    return target_cols


def _normalize_column_name(name):
    """Normalize a column name to alphanumeric lowercase for fuzzy matching."""
    return "".join(ch for ch in str(name).strip().lower() if ch.isalnum())


def resolve_mapping_columns(mapping_columns):
    """Resolve required mapping columns even with spacing/case differences.

    Example: " Updated Sample Name " and "updatedsamplename" both normalize
    to the same key and can be matched reliably.
    """
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
    """Ensure all required metrics were discovered before analysis proceeds."""
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
    """Compute standard error of the mean (SEM); 0 when only one replicate."""
    return x.std() / np.sqrt(len(x)) if len(x) > 1 else 0

def _sample_type_rank(sample_type):
    """Sort key for global plot ordering: Negative, Positive, Experimental."""
    sample_type_norm = str(sample_type).strip().lower()
    if "negative" in sample_type_norm:
        return 0
    if "positive" in sample_type_norm:
        return 1
    return 2


def _positive_control_rank(sample_name):
    """Secondary sort key for positive controls to keep expected visual order."""
    sample_name_norm = str(sample_name).strip().lower()
    if "a-his" in sample_name_norm:
        return 0
    if "k1-70" in sample_name_norm:
        return 1
    return 2


def calculate_mock_expression_threshold(df, target_cols):
    """Compute expression threshold = 2x mean(mock_His/mock_FLAG expression).

    Search order:
    1) "Sample Name"
    2) "True Sample Name" (fallback if no matches in #1)
    """
    expression_col = target_cols["mirfp_expression"]
    expression_values = pd.to_numeric(df[expression_col], errors="coerce")
    sample_names = df["Sample Name"].astype(str)

    def build_mock_mask(name_series):
        """True for names containing "mock" and at least one of "his"/"flag"."""
        return (
            name_series.str.contains("mock", case=False, na=False)
            & (
                name_series.str.contains("his", case=False, na=False)
                | name_series.str.contains("flag", case=False, na=False)
            )
        )

    mock_with_marker_mask = build_mock_mask(sample_names)
    if not mock_with_marker_mask.any() and "True Sample Name" in df.columns:
        true_sample_names = df["True Sample Name"].astype(str)
        mock_with_marker_mask = build_mock_mask(true_sample_names)

    mock_marker_values = expression_values[mock_with_marker_mask]
    mock_marker_mean = mock_marker_values.mean()
    if pd.isna(mock_marker_mean):
        raise ValueError(
            "Could not compute mock expression threshold because no valid mock_His/mock_FLAG "
            "expression values were found."
        )
    return float(mock_marker_mean) * 2.0


def calculate_percent_parent_threshold(plot_data):
    """Compute key-findings threshold = 2x best non-mock FLAG %Parent mean.

    If no qualifying non-mock FLAG sample exists, return None.
    """
    true_sample_names = plot_data["Sample Name"].astype(str)

    # Primary rule: use True Sample Name values containing FLAG but not mock.
    flag_non_mock = plot_data[
        true_sample_names.str.contains("flag", case=False, na=False)
        & ~true_sample_names.str.contains("mock", case=False, na=False)
    ]
    if not flag_non_mock.empty:
        flag_group_means = (
            flag_non_mock.groupby("Sample Name", sort=False)["percent_parent_mean"]
            .mean()
            .sort_values(ascending=False)
        )
        if not flag_group_means.empty:
            return float(flag_group_means.iloc[0]) * 2.0

    return None


def calculate_percent_parent_plot_threshold(plot_data):
    """Compute plotting threshold for %Parent horizontal reference line.

    Preference:
    1) Use non-mock FLAG threshold (same logic as key findings)
    2) Fallback to strongest mock group threshold for visual guidance
    """
    # Prefer FLAG-based threshold for plotting.
    flag_threshold = calculate_percent_parent_threshold(plot_data)
    if flag_threshold is not None:
        return flag_threshold

    # Fallback rule: when only mock-labeled controls are present, choose the mock group
    # with the higher mean and use it to derive the threshold.
    true_sample_names = plot_data["Sample Name"].astype(str)
    mock_rows = plot_data[true_sample_names.str.contains("mock", case=False, na=False)]
    if mock_rows.empty:
        return None

    mock_group_means = (
        mock_rows.groupby("Sample Name", sort=False)["percent_parent_mean"]
        .mean()
        .sort_values(ascending=False)
    )
    if mock_group_means.empty:
        return None
    return float(mock_group_means.iloc[0]) * 2.0


def generate_plots(df, target_cols, output_dir, mock_expression_threshold, export_png=False):
    """Create and save all summary figures; return plot metadata + grouped data.

    Notes:
    - `export_png` is retained for CLI compatibility; PNG export is always on.
    - Figure order and sample ordering are fixed for cross-report consistency.
    - Uses explicit matplotlib bar positions to keep bars and SEM aligned.
    """
    validate_target_columns(target_cols, df.columns)

    # Color mapping
    color_map = {
        'Negative Control': '#D3D3D3', # Light Grey
        'Positive Control': '#FFB6C1', # Light Red
        'Experimental Sample': '#ADD8E6' # Light Blue
    }
    
    # Group by Sample Name and Type for plotting averages
    # We want to keep the order from the mapping file if possible
    # We use 'True Sample Name' and 'Sample Type' from the mapping file
    # Aggregation is at sample-level mean ± SEM across replicates.
    plot_data = df.groupby(['True Sample Name', 'Sample Type'], sort=False).agg({
        target_cols['percent_parent']: ['mean', get_sem],
        target_cols['mfi_ratio']: ['mean', get_sem],
        target_cols['mfi_af488']: ['mean', get_sem],
        target_cols['mirfp_expression']: ['mean', get_sem],
    }).reset_index()
    
    # Flatten columns
    plot_data.columns = ['Sample Name', 'Sample Type', 
                         'percent_parent_mean', 'percent_parent_sem',
                         'mfi_ratio_mean', 'mfi_ratio_sem',
                         'mfi_af488_mean', 'mfi_af488_sem',
                         'mirfp_expression_mean', 'mirfp_expression_sem']
    percent_parent_threshold = calculate_percent_parent_plot_threshold(plot_data)
    
    metrics = [
        {
            "metric_id": "percent_parent",
            "title": "Binding Competent Population (%Parent)",
            "y_label": "Mean Binding Competent Population (%Parent)",
            "filename": "percent_parent_plot.png",
        },
        {
            "metric_id": "mirfp_expression",
            "title": "Expression Level of Transfected Cells MFI_AF647",
            "y_label": "Expression Level (MFI_AF647)",
            "filename": "mirfp_expression_plot.png",
        },
        {
            "metric_id": "mfi_ratio",
            "title": "Binding Competent Population Single-Cell Ratio of AF488/AF647",
            "y_label": "MFI (AF488/AF647)",
            "filename": "mfi_ratio_plot.png",
        },
        {
            "metric_id": "mfi_af488",
            "title": "Binding Competent Population MFI_AF488",
            "y_label": "MFI (AF488)",
            "filename": "mfi_af488_plot.png",
        },
    ]

    # For figures only: omit samples with "Mock" in the name.
    figure_data = plot_data[~plot_data["Sample Name"].astype(str).str.contains("mock", case=False, na=False)].copy()

    # Ordering for all figures:
    # 1) Negative controls
    # 2) Positive controls
    # 3) Experimental samples (sorted by percent_parent_mean ascending)
    experimental_order = (
        figure_data[figure_data["Sample Type"].astype(str).str.contains("experimental", case=False, na=False)]
        .sort_values("percent_parent_mean", ascending=True)["Sample Name"]
        .tolist()
    )
    experimental_rank = {name: idx for idx, name in enumerate(experimental_order)}

    figure_data["sample_type_rank"] = figure_data["Sample Type"].apply(_sample_type_rank)
    figure_data["positive_rank"] = figure_data["Sample Name"].apply(_positive_control_rank)
    figure_data["experimental_rank"] = figure_data["Sample Name"].map(experimental_rank).fillna(-1)
    figure_data = figure_data.sort_values(
        by=["sample_type_rank", "positive_rank", "experimental_rank", "Sample Name"],
        ascending=[True, True, True, True],
    ).drop(columns=["sample_type_rank", "positive_rank", "experimental_rank"])

    plot_files = []
    
    for metric_cfg in metrics:
        metric_id = metric_cfg["metric_id"]
        fig, ax = plt.subplots(figsize=(14, 7))
        x_positions = np.arange(len(figure_data))
        means = figure_data[f"{metric_id}_mean"].to_numpy()
        sems = figure_data[f"{metric_id}_sem"].to_numpy()
        bar_colors = [color_map.get(sample_type, "#B0B0B0") for sample_type in figure_data["Sample Type"]]

        # Draw bars and error bars together to keep a strict 1:1 row-to-bar mapping.
        bars = ax.bar(
            x_positions,
            means,
            yerr=sems,
            capsize=3,
            color=bar_colors,
            edgecolor="black",
            linewidth=0.5,
            zorder=3,
        )

        ax.set_title(metric_cfg["title"])
        ax.set_xticks(x_positions)
        ax.set_xticklabels(figure_data["Sample Name"], rotation=45, ha="right", fontsize=8)
        ax.set_ylabel(metric_cfg["y_label"])
        ax.set_xlabel("Sample Name")
        y_max = float(np.nanmax(means)) if len(means) else 0.0
        if metric_id == "mirfp_expression":
            y_max = max(y_max, float(mock_expression_threshold))
        ax.set_ylim(0, (y_max * 1.3) if y_max > 0 else 1.0)
        ax.set_axisbelow(True)
        # Grid is drawn behind bars to preserve bar readability.
        ax.yaxis.grid(True, linestyle="--", linewidth=0.25, color="gray", alpha=0.7, which="major")

        if metric_id == "percent_parent" and percent_parent_threshold is not None:
            ax.axhline(
                    y=percent_parent_threshold,
                    color="red",
                    linestyle="--",
                    linewidth=1.0,
                    zorder=2,
            )
        if metric_id == "mirfp_expression":
            ax.axhline(
                y=mock_expression_threshold,
                color="red",
                linestyle="--",
                linewidth=1.0,
                zorder=2,
            )

        legend_order = ["Negative Control", "Positive Control", "Experimental Sample"]
        present_types = figure_data["Sample Type"].dropna().astype(str).unique().tolist()
        legend_handles = [
            Patch(facecolor=color_map[sample_type], edgecolor="black", label=sample_type)
            for sample_type in legend_order
            if sample_type in present_types
        ]
        bars_need_hatch = []
        if metric_id == "mirfp_expression":
            # Hatch bars below expression threshold to flag insufficient expression.
            for bar, mean_value in zip(bars, means):
                if float(mean_value) < float(mock_expression_threshold):
                    bar.set_hatch("//")
                    bars_need_hatch.append(True)
                else:
                    bars_need_hatch.append(False)
            if any(bars_need_hatch):
                legend_handles.append(
                    Patch(
                        facecolor="white",
                        edgecolor="black",
                        hatch="//",
                        label="Insufficient Expression",
                    )
                )
        if legend_handles:
            ax.legend(handles=legend_handles, title="Sample Type", loc="upper right")
        fig.tight_layout()
        
        filename = os.path.join(output_dir, metric_cfg["filename"])
        # Save as 150 dpi to reduce file size.
        fig.savefig(filename, dpi=150, bbox_inches='tight')
        plot_files.append((metric_cfg["title"], metric_cfg["filename"]))
        print(f"Exported {filename}")
            
        plt.close(fig)
        
    return plot_files, plot_data, percent_parent_threshold

def _format_sample_list(sample_names):
    """Format sample-name lists for readable markdown bullet output."""
    return ", ".join(sample_names) if sample_names else "None"


def generate_report(plot_data, plot_files, target_cols, output_dir, mock_expression_threshold, percent_parent_threshold):
    """Write markdown report with key findings, summary table, and figure links."""
    output_path = os.path.join(output_dir, "experiment_summary.md")
    with open(output_path, "w") as f:
        f.write("# Flow Cytometry Analysis Summary\n\n")

        # Subset 1: experimental samples above expression threshold.
        experimental_samples = plot_data[
            plot_data["Sample Type"].astype(str).str.contains("experimental", case=False, na=False)
        ].copy()
        passed_expression = experimental_samples[
            experimental_samples["mirfp_expression_mean"] > float(mock_expression_threshold)
        ].copy()
        passed_expression_names = passed_expression["Sample Name"].astype(str).tolist()

        # Subset 2: from subset 1, samples above 2x FLAG %Parent threshold.
        key_findings_flag_threshold = calculate_percent_parent_threshold(plot_data)
        if key_findings_flag_threshold is None:
            passed_percent_parent = passed_expression.iloc[0:0].copy()
        else:
            passed_percent_parent = passed_expression[
                passed_expression["percent_parent_mean"] > float(key_findings_flag_threshold)
            ].copy()
        passed_percent_parent_names = passed_percent_parent["Sample Name"].astype(str).tolist()

        controls = plot_data[
            plot_data["Sample Type"].astype(str).str.contains("negative|positive", case=False, na=False)
        ]
        controls_ratio_mean = float(controls["mfi_ratio_mean"].mean()) if not controls.empty else float("nan")
        # Subset 3: from subset 2, samples above average control ratio.
        passed_ratio = passed_percent_parent[
            passed_percent_parent["mfi_ratio_mean"] > controls_ratio_mean
        ].copy()
        passed_ratio_names = passed_ratio["Sample Name"].astype(str).tolist()

        f.write("## Key Findings\n\n")
        f.write(
            f"- Experimental samples >2X mock expression: {len(passed_expression_names)} "
            f"({_format_sample_list(passed_expression_names)})\n"
        )
        if key_findings_flag_threshold is not None:
            f.write(
                f"- From that subset, samples >2X FLAG binding %Parent threshold: {len(passed_percent_parent_names)} "
                f"({_format_sample_list(passed_percent_parent_names)})\n"
            )
        else:
            f.write(
                "- From that subset, samples >2X FLAG binding %Parent threshold: "
                "threshold unavailable (no qualifying FLAG non-mock control found)\n"
            )
        if np.isnan(controls_ratio_mean):
            f.write(
                "- From that subset, samples above mean AF488/AF647 ratio of all controls: "
                "control average unavailable\n\n"
            )
        else:
            f.write(
                f"- From that subset, samples above mean AF488/AF647 ratio of all controls "
                f"({controls_ratio_mean:.2f}): {len(passed_ratio_names)} "
                f"({_format_sample_list(passed_ratio_names)})\n\n"
            )
        
        f.write("## Data Table\n\n")
        # Format table for markdown
        table_df = plot_data.copy()
        
        # Map back to original column names for the table headers if desired, 
        # but let's use cleaner names for the table.
        display_cols = {
            'Sample Name': 'Sample Name',
            'Sample Type': 'Sample Type',
            'mock_expression_pass': '>2X Mock Expression',
            'mirfp_expression': 'Expression Level (MFI_AF647)',
            'percent_parent': 'Singlets/AF647(+)/AF488(+) %Parent',
            'mfi_ratio': 'MFI Ratio (AF488/AF647)',
            'mfi_af488': 'MFI AF488'
        }

        table_df[display_cols['mock_expression_pass']] = table_df.apply(
            lambda row: (
                "Yes"
                if (
                    # Mock rows are intentionally marked "No" in pass/fail output.
                    "mock" not in str(row["Sample Name"]).lower()
                    and float(row["mirfp_expression_mean"]) > float(mock_expression_threshold)
                )
                else "No"
            ),
            axis=1,
        )
        
        for metric in ['mirfp_expression', 'percent_parent', 'mfi_ratio', 'mfi_af488']:
            table_df[display_cols[metric]] = table_df.apply(lambda r: f"{r[metric+'_mean']:.2f} ± {r[metric+'_sem']:.2f}", axis=1)
        
        final_table = table_df[
            ['Sample Name', 'Sample Type', display_cols['mock_expression_pass']]
            + [display_cols[m] for m in ['mirfp_expression', 'percent_parent', 'mfi_ratio', 'mfi_af488']]
        ]
        f.write(final_table.to_markdown(index=False))
        f.write("\n\n")
        
        f.write("## Figures\n\n")
        for title, png_filename in plot_files:
            f.write(f"### {title}\n")
            f.write(f"![{title}]({png_filename})\n\n")
            
    print(f"Generated {output_path}")

def main():
    """CLI entrypoint for end-to-end analysis workflow."""
    parser = argparse.ArgumentParser(description="Analyze Flow Cytometry Data")
    parser.add_argument("data_dir", help="Path to directory containing raw CSV and plate mapping CSV")
    parser.add_argument("--export-png", action="store_true", help="Export plots as PNG files")
    args = parser.parse_args()
    
    try:
        # 1) Resolve output location and clean/merge source files.
        output_dir = get_output_dir(args.data_dir)
        print(f"Writing outputs to: {output_dir}")

        merged_df = clean_and_merge(args.data_dir, output_dir)

        # 2) Resolve metric columns and compute analysis thresholds.
        target_cols = identify_columns(merged_df)
        validate_target_columns(target_cols, merged_df.columns)
        mock_expression_threshold = calculate_mock_expression_threshold(merged_df, target_cols)
        
        print("Identified target columns:")
        for k, v in target_cols.items():
            print(f"  {k}: {v}")

        # 3) Generate figures and markdown summary report.
        plot_files, plot_data, percent_parent_threshold = generate_plots(
            merged_df, target_cols, output_dir, mock_expression_threshold, args.export_png
        )
        generate_report(
            plot_data,
            plot_files,
            target_cols,
            output_dir,
            mock_expression_threshold,
            percent_parent_threshold,
        )
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
