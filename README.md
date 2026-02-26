# Flow Cytometry Analysis

This repository contains Python scripts to process FlowJo-exported CSV data, merge it with plate mapping metadata, compute experiment-level summary statistics, and generate standardized figures/reports for rapid interpretation.

## Table of Contents

- [Big-Picture Overview](#big-picture-overview)
- [Scripts](#scripts)
  - [`analyze_flow.py`](#analyze_flowpy)
  - [`analyze_flow_V2.py`](#analyze_flow_v2py)
- [Data Requirements](#data-requirements)
  - [Input Data Type](#input-data-type)
  - [Required Inputs Per Run](#required-inputs-per-run)
  - [Raw CSV Expectations](#raw-csv-expectations)
  - [Mapping CSV Expectations](#mapping-csv-expectations)
- [`analyze_flow_V2.py`: Step-by-Step Pipeline](#analyze_flow_v2py-step-by-step-pipeline)
- [Output Files Per Run](#output-files-per-run)
- [Running the Script](#running-the-script)
- [Troubleshooting](#troubleshooting)
- [Quick Diagnostic Checklist](#quick-diagnostic-checklist)

## Big-Picture Overview

The workflow is designed for experiments where each row is a sample (often with replicates), and FlowJo exports many metric columns per sample. The scripts:

1. Find the required raw CSV and plate mapping CSV.
2. Clean known artifact rows from raw data (`Mean`, `SD`).
3. Merge mapping metadata (true sample names, control/experimental labels, replicate IDs).
4. Compute per-sample means and SEMs for key metrics.
5. Compute threshold lines used for interpretation.
6. Export a cleaned intermediate CSV, a markdown summary, and 4 PNG figures.

## Scripts

### `analyze_flow.py`
- Original/primary implementation of the analysis pipeline.
- Includes robust handling for:
  - input file discovery,
  - mapping-column normalization/validation,
  - threshold calculations,
  - plotting and markdown report generation.
- Produces standardized outputs in a project-local folder named `<input_folder>_analyzed_data`.

### `analyze_flow_V2.py`
- Optimized/refactored standalone version with the same analysis behavior and output structure.
- Organized into clearer helper functions for easier debugging and future optimization:
  - file discovery,
  - merge/clean logic,
  - aggregation,
  - plotting,
  - report writing.

## Data Requirements

### Input Data Type
- CSV files exported from FlowJo and experiment plate mapping sheets.

### Required Inputs Per Run
- One **raw flow data CSV**.
- One **plate mapping CSV**.

### Raw CSV Expectations
- First column contains sample identifiers (sample names).
- Remaining columns contain flow metrics.
- Artifact rows labeled `Mean` and `SD` may be present and are removed.
- Must include columns that can be pattern-matched to these metrics:
  - `percent_parent`: contains `Freq. of Parent (%)` and `AF488(+)`
  - `mfi_ratio`: contains `Ratio_AF488_AF647`
  - `mfi_af488`: contains `AF488-A` and `AF488(+)`
  - `mirfp_expression`: contains `a-FLAG_AF647(+)` and `Geometric Mean (R1-A :: miRFP-A)` and does **not** contain `/a-His_AF488(+)`

### Mapping CSV Expectations
- Must provide columns (case/spacing tolerant):
  - `Sample Name`
  - `Updated Sample Name` (renamed to `True Sample Name` in outputs)
  - `Sample Type` (e.g., Negative Control, Positive Control, Experimental Sample)
  - `Replicate`

## `analyze_flow_V2.py`: Step-by-Step Pipeline

1. **Resolve Output Directory**
   - Creates `<input_folder>_analyzed_data` in the project root if needed.

2. **Find Input CSVs**
   - Searches the provided directory for:
     - mapping file name containing `plate_mapping`
     - raw file names containing `flowjo table`, `ratio_reanalysis`, or `reanaly`
   - If not found at top level, searches one nested directory level.

3. **Load and Clean Raw Data**
   - Reads raw CSV.
   - Treats first column as sample identifier.
   - Removes rows where sample identifier is `Mean` or `SD` (case-insensitive).
   - Drops rows with empty sample identifiers.

4. **Validate and Normalize Mapping Columns**
   - Reads mapping CSV.
   - Normalizes header names for robust matching.
   - Renames mapping columns to canonical names used downstream.

5. **Merge Raw + Mapping Data**
   - Left-joins mapping metadata onto raw rows by sample name.
   - Fails fast if any row is missing mapped metadata (`True Sample Name`, `Sample Type`, or `Replicate`).
   - Removes duplicate/unused columns (including `Unnamed:*`).
   - Front-loads metadata columns:
     - `Sample Name`, `True Sample Name`, `Sample Type`, `Replicate`
   - Writes `processed_flow_data.csv`.

6. **Detect Required Metric Columns**
   - Pattern-matches required metric columns.
   - Stops with explicit error if any required metric cannot be identified.

7. **Compute Thresholds**
   - **Expression threshold (`mock_expression_threshold`)**
     - `2 x mean` of `mirfp_expression` values for rows with names containing:
       - `mock` and (`his` or `flag`) (case-insensitive)
     - Checks `Sample Name` first; falls back to `True Sample Name`.
   - **Binding competent threshold for key findings (`percent_parent_threshold`)**
     - `2 x` highest group mean `%Parent` among non-mock rows containing `flag`.
   - **Plot-only `%Parent` fallback threshold**
     - If no qualifying non-mock `flag` group exists, uses `2 x` highest mock-group `%Parent` mean for the red line on the `%Parent` plot.

8. **Aggregate for Plot/Table**
   - Groups by `True Sample Name` + `Sample Type`.
   - Computes mean and SEM for:
     - `%Parent`
     - `MFI ratio (AF488/AF647)`
     - `MFI AF488`
     - `Expression (MFI_AF647 / miRFP metric)`

9. **Prepare Figure Ordering**
   - Excludes `mock` samples from figures.
   - Orders all figures consistently:
     - Negative controls
     - Positive controls (with specific positive-control sub-order)
     - Experimental samples sorted ascending by `%Parent` mean

10. **Generate 4 Figures**
    - `%Parent` plot
    - Expression plot (MFI_AF647)
    - Ratio plot (AF488/AF647)
    - MFI_AF488 plot
    - Styling includes:
      - fixed color mapping by sample type,
      - mean ± SEM bars,
      - dashed y-grid,
      - threshold lines,
      - hatching for expression bars below threshold.
    - Saves PNGs at `150 dpi`.

11. **Generate Markdown Report**
    - Writes `experiment_summary.md` with:
      - **Key Findings** section (3-step subset filtering),
      - data table with formatted `mean ± sem`,
      - `>2X Mock Expression` pass/fail column,
      - figure links to local PNG files (not base64).

## Output Files Per Run

Inside `<input_folder>_analyzed_data`:
- `processed_flow_data.csv`
- `experiment_summary.md`
- `percent_parent_plot.png`
- `mirfp_expression_plot.png`
- `mfi_ratio_plot.png`
- `mfi_af488_plot.png`

## Running the Script

Example:

```bash
python3 analyze_flow_V2.py "/absolute/path/to/data_folder"
```

Optional flag (retained for CLI compatibility):

```bash
python3 analyze_flow_V2.py "/absolute/path/to/data_folder" --export-png
```

Note: PNG export is already performed by default in the current implementation.

## Troubleshooting

### 1) `Could not find required CSV files in ...`
- **Cause:** The script could not locate both required CSVs (`raw` + `plate_mapping`) in the provided folder (or one nested level below).
- **Check:**
  - A mapping CSV filename includes `plate_mapping`.
  - The raw CSV filename includes one of: `flowjo table`, `ratio_reanalysis`, or `reanaly`.
- **Fix:** Point `data_dir` to the folder that directly contains those files (or adjust file naming to match expected patterns).

### 2) `Mapping CSV is missing required columns...`
- **Cause:** Mapping CSV is missing one or more required metadata fields.
- **Required fields (name matching is case/spacing tolerant):**
  - `Sample Name`
  - `Updated Sample Name`
  - `Sample Type`
  - `Replicate`
- **Fix:** Add or rename mapping columns so all four required fields are present.

### 3) `Found raw samples missing mapping metadata...`
- **Cause:** Some sample names in the raw CSV do not match mapping rows, so merged metadata is missing.
- **Check:**
  - Raw sample IDs match mapping `Sample Name` values exactly (aside from casing/formatting differences not automatically normalized).
  - No accidental extra spaces/suffixes in one file but not the other.
- **Fix:** Correct mapping rows or raw sample IDs so every raw sample maps to a `True Sample Name`, `Sample Type`, and `Replicate`.

### 4) `Unable to identify all required metric columns...`
- **Cause:** One or more required FlowJo metric columns were not found in the raw CSV.
- **Fix:** Verify raw export includes all required metric outputs and that column headers still contain expected patterns:
  - `%Parent` (`Freq. of Parent (%)` + `AF488(+)`)
  - ratio (`Ratio_AF488_AF647`)
  - AF488 MFI (`AF488-A` + `AF488(+)`)
  - expression (`a-FLAG_AF647(+)` + `Geometric Mean (R1-A :: miRFP-A)` excluding `/a-His_AF488(+)`)

### 5) `Could not compute mock expression threshold...`
- **Cause:** No rows matched mock expression controls (`mock` and `his/flag`) with valid numeric expression values.
- **Behavior:** Script checks `Sample Name` first, then `True Sample Name`.
- **Fix:** Ensure mock control naming includes `mock` plus `his` or `flag`, and that expression metric values are present/numeric.

### 6) Plots look different between runs
- **Possible reasons:**
  - Input CSVs changed (including mapping edits).
  - Different dataset path was used.
  - Existing output folder was overwritten by a newer run.
- **Fix:** Re-run with the same input folder and compare `processed_flow_data.csv` first, then markdown + plots.

### 7) Script runs but output folder is not where expected
- **Behavior:** Output is always written to the project folder (same folder as script), named `<input_folder>_analyzed_data`.
- **Fix:** Check the project root for the generated folder name derived from the basename of your `data_dir`.

## Quick Diagnostic Checklist

Use this sequence when a run fails or outputs look unexpected:

1. **Confirm path correctness**
   - Ensure `data_dir` points to the intended dataset folder.

2. **Verify required CSVs exist**
   - Check for one raw CSV and one mapping CSV (`plate_mapping` in filename).

3. **Check mapping schema**
   - Confirm mapping contains: `Sample Name`, `Updated Sample Name`, `Sample Type`, `Replicate`.

4. **Check sample-name alignment**
   - Verify raw sample IDs match mapping `Sample Name` values.

5. **Validate metric columns in raw export**
   - Confirm required metric patterns are present (`%Parent`, ratio, AF488 MFI, expression metric).

6. **Run once and inspect first failing message**
   - Most failures are explicit and point directly to missing files, columns, or mappings.

7. **Inspect generated `processed_flow_data.csv` first**
   - If outputs look wrong, validate merged metadata + raw values before reviewing plots.

8. **Re-run on same path to confirm reproducibility**
   - Output folder is overwritten; compare markdown/table values across runs.
