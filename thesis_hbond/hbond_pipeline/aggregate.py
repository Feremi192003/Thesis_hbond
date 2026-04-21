from __future__ import annotations

"""
Helpers for averaging and reshaping pipeline outputs.

This module takes file-level, residue-level metric tables and conversts them into condition/isoform summaries used in other modules in the pipeline.
"""

import pandas as pd


def aggregate_metrics(file_metrics_df: pd.DataFrame) -> pd.DataFrame:
    """Average file-level metrics across replicates for each condition and isoform pair.

    Takes table containing one row per processed CSV file, may include Core metrics + residue-specific columns.

    Returns:
        pd.DataFrame, a table with one row per condition and isoform pair and mean values for every numeric metric column of interest.
    """
    #Start with the main-file level metrics, append any residue count columns that were created dynamically during metric computatiion.
    numeric_cols = ["sis", "nucleotide_sis", "trt", "ehbs"] + [
        c for c in file_metrics_df.columns if c.startswith("count_")
    ]

    #Group replicate rows by condition, isoform, average the metrics, then sort so tables and plots have stable ordering.
    aggregated_df = (
        file_metrics_df.groupby(["condition", "isoform"], as_index=False)[numeric_cols]
        .mean()
        .sort_values(["condition", "isoform"])
        .reset_index(drop=True)
    )

    return aggregated_df

def aggregate_residue_metrics(file_residue_metrics_df: pd.DataFrame) -> pd.DataFrame:
    """Average residue-level contributions across replicates.

   Takes a table with one row per residue file with metrics.

    Returns:
        pd.DataFrame: which is a summary indexed by wt/y220c and iso, if input is empty, then empty table returned so the rest of pipeline still runs.
    """
    #preserve order.
    if file_residue_metrics_df.empty:
        return pd.DataFrame(
            columns=[
                "condition",
                "isoform",
                "residue_label",
                "canonical_specific_residue",
                "canonical_nucleotide_residue",
                "is_specific_residue",
                "is_nucleotide_residue",
                "total_occupancy",
                "unique_dna_contacts",
                "sis_contribution",
                "nucleotide_sis_contribution",
            ]
        )

    #numeric res quantities that should be averaged.
    numeric_cols = [
        "total_occupancy",
        "unique_dna_contacts",
        "sis_contribution",
        "nucleotide_sis_contribution",
    ]

    #Identify a residue within a given condition and isoform context that must be preserved while aggregating.
    group_cols = [
        "condition",
        "isoform",
        "residue_label",
        "canonical_specific_residue",
        "canonical_nucleotide_residue",
        "is_specific_residue",
        "is_nucleotide_residue",
    ]

    aggregated_residue_df = (
        file_residue_metrics_df.groupby(group_cols, as_index=False)[numeric_cols]
        .mean()
        .sort_values(["condition", "isoform", "residue_label"])
        .reset_index(drop=True)
    )

    return aggregated_residue_df


def condition_comparison_table(aggregated_df: pd.DataFrame) -> pd.DataFrame:
    """Build a wt/y220c summary table.

    collapses the table into a single row per contribution and extracts SIS, nucSis, TRT, EHBS to compare for both systems.

    """
    #Sum across isoforms so each condition is a single total 
    summary = (
        aggregated_df.groupby("condition", as_index=False)[["sis", "nucleotide_sis", "trt", "ehbs"]]
        .sum()
    )

    wt_row = summary.loc[summary["condition"] == "WT"]
    mut_row = summary.loc[summary["condition"] == "Y220C"]
    
    #Using 0.0 fallbacks so table can still be made even if one condition is missing.
    wt_sis = float(wt_row["sis"].iloc[0]) if not wt_row.empty else 0.0
    wt_nucleotide_sis = float(wt_row["nucleotide_sis"].iloc[0]) if not wt_row.empty else 0.0
    wt_trt = float(wt_row["trt"].iloc[0]) if not wt_row.empty else 0.0
    wt_ehbs = float(wt_row["ehbs"].iloc[0]) if not wt_row.empty else 0.0

    mut_sis = float(mut_row["sis"].iloc[0]) if not mut_row.empty else 0.0
    mut_nucleotide_sis = float(mut_row["nucleotide_sis"].iloc[0]) if not mut_row.empty else 0.0
    mut_trt = float(mut_row["trt"].iloc[0]) if not mut_row.empty else 0.0
    mut_ehbs = float(mut_row["ehbs"].iloc[0]) if not mut_row.empty else 0.0

    comparison_df = pd.DataFrame(
        {
            "measure": [
                "Specific Interaction Score",
                "Nucleotide Specific Interaction Score",
                "Total Residence Time",
                "Exhaustive Hydrogen Bond Sum",
            ],
            "WT": [wt_sis, wt_nucleotide_sis, wt_trt, wt_ehbs],
            "Y220C": [mut_sis, mut_nucleotide_sis, mut_trt, mut_ehbs],
        }
    )

    return comparison_df


def build_condition_tables(aggregated_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split summary into separated wt, y220c tables.

    takes isoform summary created by aggregate_metrics

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: two tables for isoforms, one for wt and one for y220c.
    """
    #keep only columns that are exported for condition-specific summaries and sort by isoform for easier reading.
    wt_df = (
        aggregated_df.loc[
            aggregated_df["condition"] == "WT",
            ["isoform", "sis", "nucleotide_sis", "trt", "ehbs"],
        ]
        .sort_values("isoform")
        .reset_index(drop=True)
    )

    y220c_df = (
        aggregated_df.loc[
            aggregated_df["condition"] == "Y220C",
            ["isoform", "sis", "nucleotide_sis", "trt", "ehbs"],
        ]
        .sort_values("isoform")
        .reset_index(drop=True)
    )

    return wt_df, y220c_df