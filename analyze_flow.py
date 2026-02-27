"""Optimized standalone V2 flow cytometry analysis script.

This version keeps behavior/output equivalent to `analyze_flow.py` while
reducing duplication and improving maintainability.

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

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

REQUIRED_METRIC_KEYS = ("percent_parent", "mfi_ratio", "mfi_af488", "mirfp_expression")
REQUIRED_MAPPING_COLUMNS = ("Sample Name", "Updated Sample Name", "Sample Type", "Replicate")

COLOR_MAP = {
    "Negative Control": "#D3D3D3",
    "Positive Control": "#FFB6C1",
    "Experimental Sample": "#ADD8E6",
}

METRIC_CONFIGS = [
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


def get_output_dir(data_dir):
    """Build the project-local output directory for a given input dataset folder."""
    input_dir = os.path.abspath(os.path.normpath(data_dir))
    input_name = os.path.basename(input_dir)
    project_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(project_dir, f"{input_name}_analyzed_data")
    os.makedirs(output_dir, exist_ok=True)
    return output_dir


def _pick_csvs_from_dir(dir_path):
    """Find candidate raw + mapping CSV files in a single directory."""
    raw_csv = None
    mapping_csv = None
    for name in os.listdir(dir_path):
        if not name.endswith(".csv"):
            continue
        lower_name = name.lower()
        if "plate_mapping" in lower_name and mapping_csv is None:
            mapping_csv = os.path.join(dir_path, name)
        elif raw_csv is None and (
            "flowjo table" in lower_name
            or "ratio_reanalysis" in lower_name
            or "reanaly" in lower_name
        ):
            raw_csv = os.path.join(dir_path, name)
    return raw_csv, mapping_csv


def find_input_csvs(data_dir):
    """Locate required CSVs in `data_dir` or one nested level below it."""
    # First, try the provided directory directly.
    raw_csv, mapping_csv = _pick_csvs_from_dir(data_dir)
    if raw_csv and mapping_csv:
        return raw_csv, mapping_csv

    # Fallback: search immediate child folders (common for archived exports).
    for child in os.listdir(data_dir):
        child_path = os.path.join(data_dir, child)
        if not os.path.isdir(child_path):
            continue
        raw_csv, mapping_csv = _pick_csvs_from_dir(child_path)
        if raw_csv and mapping_csv:
            return raw_csv, mapping_csv

    raise FileNotFoundError(f"Could not find required CSV files in {data_dir}")


def _normalize_column_name(name):
    """Normalize column names so minor formatting differences still match."""
    return "".join(ch for ch in str(name).strip().lower() if ch.isalnum())


def resolve_mapping_columns(mapping_columns):
    """Resolve mapping headers to canonical names expected by downstream logic."""
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

    # Build an actionable error message when required mapping metadata is missing.
    if missing:
        available_preview = ", ".join([str(col) for col in mapping_columns])
        raise ValueError(
            "Mapping CSV is missing required columns. "
            f"Missing: {', '.join(missing)}. "
            f"Available columns: {available_preview}"
        )

    return resolved


def clean_and_merge(data_dir, output_dir):
    """Load raw/mapping CSVs, clean artifacts, merge metadata, and export merged CSV."""
    raw_csv, mapping_csv = find_input_csvs(data_dir)
    print(f"Loading raw data from: {raw_csv}")
    print(f"Loading mapping from: {mapping_csv}")

    # The first raw column is treated as sample identifier regardless of its header text.
    raw_df = pd.read_csv(raw_csv)
    sample_col = raw_df.columns[0]

    # Remove FlowJo artifact summary rows and empty sample rows.
    normalized_sample_names = raw_df[sample_col].astype(str).str.strip().str.lower()
    raw_df = raw_df[~normalized_sample_names.isin(["mean", "sd"])]
    raw_df = raw_df.dropna(subset=[sample_col])

    # Standardize mapping column names so all downstream code uses stable keys.
    mapping_df = pd.read_csv(mapping_csv)
    mapping_col_map = resolve_mapping_columns(mapping_df.columns)
    mapping_df = mapping_df.rename(
        columns={
            mapping_col_map["Sample Name"]: "Sample Name",
            mapping_col_map["Updated Sample Name"]: "Updated Sample Name",
            mapping_col_map["Sample Type"]: "Sample Type",
            mapping_col_map["Replicate"]: "Replicate",
        }
    )

    # Left join keeps every raw sample row; we fail fast on missing mapping metadata below.
    merged_df = pd.merge(
        raw_df,
        mapping_df,
        left_on=sample_col,
        right_on="Sample Name",
        how="left",
    )

    # Explicitly fail on unmatched samples to avoid silent drops in later groupby operations.
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

    # Keep the raw sample identifier column first and avoid duplicate merge-key columns.
    if sample_col != "Sample Name":
        merged_df = merged_df.drop(columns=["Sample Name"], errors="ignore")
        merged_df = merged_df.rename(columns={sample_col: "Sample Name"})

    # User-facing naming convention.
    merged_df = merged_df.rename(columns={"Updated Sample Name": "True Sample Name"})

    # Drop FlowJo placeholder columns with no analysis value.
    drop_unnamed_cols = [col for col in merged_df.columns if str(col).startswith("Unnamed:")]
    if drop_unnamed_cols:
        merged_df = merged_df.drop(columns=drop_unnamed_cols, errors="ignore")

    # Front-load metadata columns to make exported CSV easier to inspect manually.
    leading_cols = ["Sample Name", "True Sample Name", "Sample Type", "Replicate"]
    remaining_cols = [col for col in merged_df.columns if col not in leading_cols]
    merged_df = merged_df[leading_cols + remaining_cols]

    output_path = os.path.join(output_dir, "processed_flow_data.csv")
    merged_df.to_csv(output_path, index=False)
    print(f"Saved merged data to {output_path}")

    return merged_df


def identify_columns(df):
    """Discover required metric columns using robust pattern matching."""
    cols = df.columns
    target_cols = {}

    parent_match = [c for c in cols if "Freq. of Parent (%)" in str(c) and "AF488(+)" in str(c)]
    ratio_match = [c for c in cols if "Ratio_AF488_AF647" in str(c)]
    af488_match = [c for c in cols if "AF488-A" in str(c) and "AF488(+)" in str(c)]
    mirfp_match = [
        c
        for c in cols
        if "a-FLAG_AF647(+)" in str(c)
        and "Geometric Mean (R1-A :: miRFP-A)" in str(c)
        and "/a-His_AF488(+)" not in str(c)
    ]

    # Use the last match to prefer the most specific column when multiple match patterns.
    if parent_match:
        target_cols["percent_parent"] = parent_match[-1]
    if ratio_match:
        target_cols["mfi_ratio"] = ratio_match[-1]
    if af488_match:
        target_cols["mfi_af488"] = af488_match[-1]
    if mirfp_match:
        target_cols["mirfp_expression"] = mirfp_match[-1]

    return target_cols


def validate_target_columns(target_cols, available_columns):
    """Verify all required metric keys were discovered before analysis proceeds."""
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
    """Compute standard error of the mean; return 0 for singleton groups."""
    return x.std() / np.sqrt(len(x)) if len(x) > 1 else 0


def _sample_type_rank(sample_type):
    """Sort rank to keep plot groups ordered: Negative -> Positive -> Experimental."""
    sample_type_norm = str(sample_type).strip().lower()
    if "negative" in sample_type_norm:
        return 0
    if "positive" in sample_type_norm:
        return 1
    return 2


def _positive_control_rank(sample_name):
    """Secondary ordering within positives to keep controls visually consistent."""
    sample_name_norm = str(sample_name).strip().lower()
    if "a-his" in sample_name_norm:
        return 0
    if "k1-70" in sample_name_norm:
        return 1
    return 2


def calculate_mock_expression_threshold(df, target_cols):
    """Compute expression threshold: 2x mean of mock_His/mock_FLAG expression."""
    expression_col = target_cols["mirfp_expression"]
    expression_values = pd.to_numeric(df[expression_col], errors="coerce")

    def build_mock_mask(name_series):
        # Case-insensitive marker filter requested by user.
        return (
            name_series.str.contains("mock", case=False, na=False)
            & (
                name_series.str.contains("his", case=False, na=False)
                | name_series.str.contains("flag", case=False, na=False)
            )
        )

    # Fallback to True Sample Name if mock markers are absent in raw Sample Name.
    mock_with_marker_mask = build_mock_mask(df["Sample Name"].astype(str))
    if not mock_with_marker_mask.any() and "True Sample Name" in df.columns:
        mock_with_marker_mask = build_mock_mask(df["True Sample Name"].astype(str))

    mock_marker_mean = expression_values[mock_with_marker_mask].mean()
    if pd.isna(mock_marker_mean):
        raise ValueError(
            "Could not compute mock expression threshold because no valid mock_His/mock_FLAG "
            "expression values were found."
        )
    return float(mock_marker_mean) * 2.0


def calculate_percent_parent_threshold(plot_data):
    """Compute 2x threshold from strongest non-mock FLAG %Parent control group."""
    sample_names = plot_data["Sample Name"].astype(str)
    flag_non_mock = plot_data[
        sample_names.str.contains("flag", case=False, na=False)
        & ~sample_names.str.contains("mock", case=False, na=False)
    ]
    if flag_non_mock.empty:
        return None

    flag_group_means = (
        flag_non_mock.groupby("Sample Name", sort=False)["percent_parent_mean"]
        .mean()
        .sort_values(ascending=False)
    )
    if flag_group_means.empty:
        return None
    return float(flag_group_means.iloc[0]) * 2.0


def calculate_percent_parent_plot_threshold(plot_data):
    """Plot threshold for %Parent; fallback to strongest mock group when needed."""
    flag_threshold = calculate_percent_parent_threshold(plot_data)
    if flag_threshold is not None:
        return flag_threshold

    sample_names = plot_data["Sample Name"].astype(str)
    mock_rows = plot_data[sample_names.str.contains("mock", case=False, na=False)]
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


def build_plot_data(df, target_cols):
    """Aggregate replicate-level rows into sample-level means and SEMs."""
    grouped = (
        df.groupby(["True Sample Name", "Sample Type"], sort=False)
        .agg(
            {
                target_cols["percent_parent"]: ["mean", get_sem],
                target_cols["mfi_ratio"]: ["mean", get_sem],
                target_cols["mfi_af488"]: ["mean", get_sem],
                target_cols["mirfp_expression"]: ["mean", get_sem],
            }
        )
        .reset_index()
    )
    # Flatten MultiIndex columns created by grouped aggregation.
    grouped.columns = [
        "Sample Name",
        "Sample Type",
        "percent_parent_mean",
        "percent_parent_sem",
        "mfi_ratio_mean",
        "mfi_ratio_sem",
        "mfi_af488_mean",
        "mfi_af488_sem",
        "mirfp_expression_mean",
        "mirfp_expression_sem",
    ]
    return grouped


def build_figure_data(plot_data):
    """Filter/sort rows used for all figures to keep ordering consistent."""
    # Figures omit mock rows by requirement, while tables keep them.
    figure_data = plot_data[
        ~plot_data["Sample Name"].astype(str).str.contains("mock", case=False, na=False)
    ].copy()

    # Experimental bars are globally ordered by percent_parent_mean and reused across all plots.
    experimental_order = (
        figure_data[
            figure_data["Sample Type"].astype(str).str.contains(
                "experimental", case=False, na=False
            )
        ]
        .sort_values("percent_parent_mean", ascending=True)["Sample Name"]
        .tolist()
    )
    experimental_rank = {name: idx for idx, name in enumerate(experimental_order)}

    figure_data["sample_type_rank"] = figure_data["Sample Type"].apply(_sample_type_rank)
    figure_data["positive_rank"] = figure_data["Sample Name"].apply(_positive_control_rank)
    figure_data["experimental_rank"] = (
        figure_data["Sample Name"].map(experimental_rank).fillna(-1)
    )

    return figure_data.sort_values(
        by=["sample_type_rank", "positive_rank", "experimental_rank", "Sample Name"],
        ascending=[True, True, True, True],
    ).drop(columns=["sample_type_rank", "positive_rank", "experimental_rank"])


def draw_metric_plot(
    ax,
    figure_data,
    metric_cfg,
    mock_expression_threshold,
    percent_parent_threshold,
):
    """Render one bar chart (mean ± SEM) with requested styling and thresholds."""
    metric_id = metric_cfg["metric_id"]
    x_positions = np.arange(len(figure_data))
    means = figure_data[f"{metric_id}_mean"].to_numpy()
    sems = figure_data[f"{metric_id}_sem"].to_numpy()
    bar_colors = [COLOR_MAP.get(st, "#B0B0B0") for st in figure_data["Sample Type"]]

    # Explicit bar + yerr keeps data-to-bar mapping deterministic.
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

    # Scale axis with headroom; ensure threshold line fits on mirfp expression plot.
    y_max = float(np.nanmax(means)) if len(means) else 0.0
    if metric_id == "mirfp_expression":
        y_max = max(y_max, float(mock_expression_threshold))
    ax.set_ylim(0, (y_max * 1.3) if y_max > 0 else 1.0)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, linestyle="--", linewidth=0.25, color="gray", alpha=0.7, which="major")

    # Draw requested threshold reference lines.
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

    # Compose legend only for sample types present in the filtered figure data.
    legend_order = ["Negative Control", "Positive Control", "Experimental Sample"]
    present_types = figure_data["Sample Type"].dropna().astype(str).unique().tolist()
    legend_handles = [
        Patch(facecolor=COLOR_MAP[sample_type], edgecolor="black", label=sample_type)
        for sample_type in legend_order
        if sample_type in present_types
    ]

    # Hatch insufficient-expression bars and add legend entry when applicable.
    if metric_id == "mirfp_expression":
        has_hatched = False
        for bar, mean_value in zip(bars, means):
            if float(mean_value) < float(mock_expression_threshold):
                bar.set_hatch("//")
                has_hatched = True
        if has_hatched:
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


def generate_plots(df, target_cols, output_dir, mock_expression_threshold, export_png=False):
    """Generate and save all figures; return filenames + aggregated plot dataset."""
    # Retained for CLI compatibility; this script always exports PNGs.
    del export_png
    validate_target_columns(target_cols, df.columns)

    # Shared sample ordering across all metrics prevents visual reordering confusion.
    plot_data = build_plot_data(df, target_cols)
    figure_data = build_figure_data(plot_data)
    percent_parent_threshold = calculate_percent_parent_plot_threshold(plot_data)

    # Iterate metric configuration so title/axis/file naming stays centralized.
    plot_files = []
    for metric_cfg in METRIC_CONFIGS:
        fig, ax = plt.subplots(figsize=(14, 7))
        draw_metric_plot(
            ax,
            figure_data,
            metric_cfg,
            mock_expression_threshold,
            percent_parent_threshold,
        )
        fig.tight_layout()

        filename = os.path.join(output_dir, metric_cfg["filename"])
        fig.savefig(filename, dpi=150, bbox_inches="tight")
        print(f"Exported {filename}")
        plt.close(fig)

        plot_files.append((metric_cfg["title"], metric_cfg["filename"]))

    return plot_files, plot_data, percent_parent_threshold


def _format_sample_list(sample_names):
    """Human-readable sample list for markdown bullet points."""
    return ", ".join(sample_names) if sample_names else "None"


def generate_report(
    plot_data,
    plot_files,
    target_cols,
    output_dir,
    mock_expression_threshold,
    percent_parent_threshold,
):
    """Build markdown report: key findings, summary table, and figure references."""
    # Kept in signature for compatibility with caller shape.
    del target_cols, percent_parent_threshold

    output_path = os.path.join(output_dir, "experiment_summary.md")
    with open(output_path, "w") as f:
        f.write("# Flow Cytometry Analysis Summary\n\n")

        # Key Findings subset 1: experimental samples above expression threshold.
        experimental_samples = plot_data[
            plot_data["Sample Type"].astype(str).str.contains("experimental", case=False, na=False)
        ].copy()
        passed_expression = experimental_samples[
            experimental_samples["mirfp_expression_mean"] > float(mock_expression_threshold)
        ].copy()
        passed_expression_names = passed_expression["Sample Name"].astype(str).tolist()

        # Key Findings subset 2: from subset 1, samples above 2x FLAG %Parent threshold.
        key_findings_flag_threshold = calculate_percent_parent_threshold(plot_data)
        if key_findings_flag_threshold is None:
            passed_percent_parent = passed_expression.iloc[0:0].copy()
        else:
            passed_percent_parent = passed_expression[
                passed_expression["percent_parent_mean"] > float(key_findings_flag_threshold)
            ].copy()
        passed_percent_parent_names = passed_percent_parent["Sample Name"].astype(str).tolist()

        # Key Findings subset 3: from subset 2, samples above average control ratio.
        controls = plot_data[
            plot_data["Sample Type"].astype(str).str.contains("negative|positive", case=False, na=False)
        ]
        controls_ratio_mean = float(controls["mfi_ratio_mean"].mean()) if not controls.empty else float("nan")
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

        # Build table columns with user-requested naming/ordering.
        f.write("## Data Table\n\n")
        table_df = plot_data.copy()
        display_cols = {
            "Sample Name": "Sample Name",
            "Sample Type": "Sample Type",
            "mock_expression_pass": ">2X Mock Expression",
            "mirfp_expression": "Expression Level (MFI_AF647)",
            "percent_parent": "Singlets/AF647(+)/AF488(+) %Parent",
            "mfi_ratio": "MFI Ratio (AF488/AF647)",
            "mfi_af488": "MFI AF488",
        }

        # Mark mock rows as "No" by design for pass/fail display.
        table_df[display_cols["mock_expression_pass"]] = table_df.apply(
            lambda row: (
                "Yes"
                if (
                    "mock" not in str(row["Sample Name"]).lower()
                    and float(row["mirfp_expression_mean"]) > float(mock_expression_threshold)
                )
                else "No"
            ),
            axis=1,
        )

        # Format aggregated values as "mean ± sem" strings for readability.
        for metric in ("mirfp_expression", "percent_parent", "mfi_ratio", "mfi_af488"):
            table_df[display_cols[metric]] = table_df.apply(
                lambda r: f"{r[metric + '_mean']:.2f} ± {r[metric + '_sem']:.2f}",
                axis=1,
            )

        final_table = table_df[
            ["Sample Name", "Sample Type", display_cols["mock_expression_pass"]]
            + [display_cols[m] for m in ("mirfp_expression", "percent_parent", "mfi_ratio", "mfi_af488")]
        ]
        f.write(final_table.to_markdown(index=False))
        f.write("\n\n")

        # Use markdown links to local PNG files (no base64 embedding).
        f.write("## Figures\n\n")
        for title, png_filename in plot_files:
            f.write(f"### {title}\n")
            f.write(f"![{title}]({png_filename})\n\n")

    print(f"Generated {output_path}")


def main():
    """CLI entrypoint for the end-to-end analysis/report generation workflow."""
    parser = argparse.ArgumentParser(description="Analyze Flow Cytometry Data")
    parser.add_argument("data_dir", help="Path to directory containing raw CSV and plate mapping CSV")
    parser.add_argument("--export-png", action="store_true", help="Export plots as PNG files")
    args = parser.parse_args()

    try:
        # 1) Resolve output location and merge cleaned data.
        output_dir = get_output_dir(args.data_dir)
        print(f"Writing outputs to: {output_dir}")

        merged_df = clean_and_merge(args.data_dir, output_dir)
        # 2) Identify metric columns and compute threshold(s).
        target_cols = identify_columns(merged_df)
        validate_target_columns(target_cols, merged_df.columns)
        mock_expression_threshold = calculate_mock_expression_threshold(merged_df, target_cols)

        print("Identified target columns:")
        for k, v in target_cols.items():
            print(f"  {k}: {v}")

        # 3) Generate figures and markdown report from aggregated data.
        plot_files, plot_data, percent_parent_threshold = generate_plots(
            merged_df,
            target_cols,
            output_dir,
            mock_expression_threshold,
            args.export_png,
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
        # Preserve full traceback for faster debugging in local runs.
        print(f"Error: {e}")
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
