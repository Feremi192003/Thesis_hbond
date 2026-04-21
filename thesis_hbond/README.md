# HBond Python Pipeline

This is a cleaned parallelized Python version of the Sheets-based hydrogen-bond workflow.

## What it computes

For each `*_DNA_RES_INTERACT.csv` file:

- **EHBS** = number of rows in the CSV
- **TRT** = sum of the `frac` column
- **SIS** = row-unique count of rows where at least one target residue appears in columns A/B/C
- **Per-residue contribution** = row-unique count for each target residue across columns A/B/C

A row only contributes **once** to SIS, even if the same residue appears in both protein columns or multiple target residues appear in the same row.

## Canonical target residues

- `SER_241`
- `ARG_283`
- `CYS_277`
- `ALA_276`
- `ARG_280`
- `ARG_273`
- `LYS_120`
- `ARG_248`

## Isoform offsets

- isoforms 1–3: `0`
- isoforms 4–6: `39`
- isoforms 7–9: `132`
- isoforms 10–12: `159`

These offsets are used to remap canonical residues to isoform-local numbering.

## Expected filename format

The parser expects filenames that begin with `WT` or `Y220C` and contain an isoform tag like `i1`, `i7`, `i12`.

Examples:

- `WT2_i1_DNA_RES_INTERACT.csv`
- `Y220C1_i12_DNA_RES_INTERACT.csv`
- `WT_r2_i4_DNA_RES_INTERACT.csv`

Replicate parsing is optional. The script looks for tags like `r1`, `r2`, `r3` if present.

## Installation

```bash
pip install pandas matplotlib
```

## Run

```bash
python -m hbond_pipeline.main /path/to/csv_directory -o /path/to/results
```

Skip plots:

```bash
python -m hbond_pipeline.main /path/to/csv_directory -o /path/to/results --no-plots
```

## Output files

- `file_level_metrics.csv` — one row per input CSV
- `isoform_condition_averages.csv` — average metrics by condition and isoform
- `wt_y220c_comparison.csv` — WT vs Y220C comparison table with deltas
- `row_annotations.csv` — one row per original interaction row, showing whether it matched target residues
- PNG plots for SIS, EHBS, TRT, and each target residue contribution

