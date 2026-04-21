from __future__ import annotations
"""Plotting for file-level and res level summaries.
    this turns tables into bar plots, heatmaps and residue charts.
    """
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def _build_mean_std_tables(
    file_metrics_df: pd.DataFrame,
    metric: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute per iso means and stand. dev for one metric, 
    returns two tables with isos on rows and conditions on cols to be passed to pandas/matplotlib.
    """
    grouped = (
        file_metrics_df.groupby(["isoform", "condition"])[metric]
        .agg(["mean", "std"])
        .reset_index()
    )

    mean_df = (
        grouped.pivot(index="isoform", columns="condition", values="mean")
        .sort_index()
        .fillna(0)
    )
    std_df = (
        grouped.pivot(index="isoform", columns="condition", values="std")
        .sort_index()
        .fillna(0)
    )

    return mean_df, std_df


def _build_clipped_yerr(
    mean_df: pd.DataFrame,
    std_df: pd.DataFrame,
) -> list:
    """
    Build asymmetric error bars so the lower bar never drops below zero.

    For nonnegative quantities, matplotlib can otherwise draw mean - std below 0
    when the standard deviation is larger than the mean.
    """
    mean_aligned, std_aligned = mean_df.align(std_df, join="left", axis=None, fill_value=0)

    #Clip lower error to the mean itself so there are never error bars below zero.
    lower_err = std_aligned.clip(upper=mean_aligned)
    upper_err = std_aligned

    return [lower_err.values.T, upper_err.values.T]


def plot_metric_by_isoform(
    aggregated_df: pd.DataFrame,
    file_metrics_df: pd.DataFrame,
    metric: str,
    output_dir: str | Path,
) -> Path:
    """Create grouped bar chart of a file-level metric across isos.
    """
    output_path = Path(output_dir) / f"{metric}_by_isoform.png"

    mean_df, std_df = _build_mean_std_tables(file_metrics_df, metric)
    yerr = _build_clipped_yerr(mean_df, std_df)

    ax = mean_df.plot(
        kind="bar",
        figsize=(10, 6),
        yerr=yerr,
        capsize=4,
    )

    ax.set_title(f"{metric.upper()} by isoform")
    ax.set_xlabel("Isoform")
    ax.set_ylabel(metric.upper())
    ax.set_ylim(bottom=0)
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()

    return output_path


def plot_residue_counts(
    aggregated_df: pd.DataFrame,
    file_metrics_df: pd.DataFrame,
    output_dir: str | Path,
) -> list[Path]:
    """Create one bar plot per residue count column.
    """
    output_dir = Path(output_dir)
    residue_columns = [col for col in aggregated_df.columns if col.startswith("count_")]

    output_paths: list[Path] = []

    for residue_column in residue_columns:
        output_path = output_dir / f"{residue_column}_by_isoform.png"

        mean_df, std_df = _build_mean_std_tables(file_metrics_df, residue_column)
        yerr = _build_clipped_yerr(mean_df, std_df)

        ax = mean_df.plot(
            kind="bar",
            figsize=(10, 6),
            yerr=yerr,
            capsize=4,
        )

        ax.set_title(f"{residue_column.replace('count_', '')} row-unique contribution")
        ax.set_xlabel("Isoform")
        ax.set_ylabel("Count")
        ax.set_ylim(bottom=0)
        plt.xticks(rotation=0)
        plt.tight_layout()
        plt.savefig(output_path, dpi=200)
        plt.close()

        output_paths.append(output_path)

    return output_paths


def plot_residue_metric_heatmap(
    aggregated_residue_df: pd.DataFrame,
    metric: str,
    output_dir: str | Path,
    residue_subset: str = "specific",
) -> Path:
    """ Create heatmaps for residue metrics."""
    output_dir = Path(output_dir)

    if aggregated_residue_df.empty:
        raise ValueError("aggregated_residue_df is empty; cannot generate heatmap.")

    if residue_subset == "specific":
        df = aggregated_residue_df.loc[aggregated_residue_df["is_specific_residue"]].copy()
        suffix = "specific"
    elif residue_subset == "nucleotide":
        df = aggregated_residue_df.loc[aggregated_residue_df["is_nucleotide_residue"]].copy()
        suffix = "nucleotide"
    else:
        raise ValueError("residue_subset must be 'specific' or 'nucleotide'")

    pivot = df.pivot_table(
        index=["residue_label"],
        columns=["isoform", "condition"],
        values=metric,
        aggfunc="mean",
        fill_value=0,
    )

    isoforms = sorted(df["isoform"].unique())
    residues = sorted(df["residue_label"].unique())

    heatmap_df = pd.DataFrame(index=residues, columns=isoforms, dtype=float)

    #fill heatmap cells with wt-y220c for that residue and iso.
    for residue_label in residues:
        for isoform in isoforms:
            wt = 0.0
            mut = 0.0

            if (isoform, "WT") in pivot.columns and residue_label in pivot.index:
                wt = pivot.loc[residue_label, (isoform, "WT")]
            if (isoform, "Y220C") in pivot.columns and residue_label in pivot.index:
                mut = pivot.loc[residue_label, (isoform, "Y220C")]

            heatmap_df.loc[residue_label, isoform] = wt - mut

    output_path = output_dir / f"{metric}_{suffix}_wt_minus_y220c_heatmap.png"

    fig, ax = plt.subplots(figsize=(12, max(4, 0.5 * len(heatmap_df))))
    im = ax.imshow(heatmap_df.values, aspect="auto")

    ax.set_title(f"WT - Y220C: {metric.replace('_', ' ')} ({suffix})")
    ax.set_xlabel("Isoform")
    ax.set_ylabel("Residue")
    ax.set_xticks(range(len(heatmap_df.columns)))
    ax.set_xticklabels(heatmap_df.columns)
    ax.set_yticks(range(len(heatmap_df.index)))
    ax.set_yticklabels(heatmap_df.index)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("WT - Y220C")

    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()

    return output_path


def plot_residue_metric_by_isoform(
    aggregated_residue_df: pd.DataFrame,
    metric: str,
    output_dir: str | Path,
    residue_subset: str = "nucleotide",
) -> Path:
    """ Create bar plots for one residue metric.
        Separate subplots can be created for each residue in chosen subset.
    """
    output_dir = Path(output_dir)

    if aggregated_residue_df.empty:
        raise ValueError("aggregated_residue_df is empty; cannot generate residue plot.")

    if residue_subset == "nucleotide":
        df = aggregated_residue_df.loc[aggregated_residue_df["is_nucleotide_residue"]].copy()
        suffix = "nucleotide"
    else:
        df = aggregated_residue_df.loc[aggregated_residue_df["is_specific_residue"]].copy()
        suffix = "specific"

    output_path = output_dir / f"{metric}_{suffix}_by_isoform.png"

    residues = sorted(df["residue_label"].unique())

    fig, axes = plt.subplots(len(residues), 1, figsize=(10, 4 * len(residues)))
    if len(residues) == 1:
        axes = [axes]

    for ax, residue in zip(axes, residues):
        sub = df[df["residue_label"] == residue]
        plot_df = (
            sub.pivot(index="isoform", columns="condition", values=metric)
            .sort_index()
            .fillna(0)
        )

        plot_df.plot(kind="bar", ax=ax)
        ax.set_title(residue)
        ax.set_xlabel("Isoform")
        ax.set_ylabel(metric.replace("_", " ").title())
        ax.set_ylim(bottom=0)
        ax.tick_params(axis="x", rotation=0)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()

    return output_path


def plot_top_residue_metric(
    aggregated_residue_df: pd.DataFrame,
    metric: str,
    output_dir: str | Path,
    condition: str = "WT",
    residue_subset: str = "specific",
    top_n: int = 10,
) -> Path:
    """ Plot top residues for one condition based on user chosen metric. """
    output_dir = Path(output_dir)

    if aggregated_residue_df.empty:
        raise ValueError("aggregated_residue_df is empty; cannot generate top-residue plot.")

    if residue_subset == "specific":
        df = aggregated_residue_df.loc[aggregated_residue_df["is_specific_residue"]].copy()
        suffix = "specific"
    elif residue_subset == "nucleotide":
        df = aggregated_residue_df.loc[aggregated_residue_df["is_nucleotide_residue"]].copy()
        suffix = "nucleotide"
    else:
        raise ValueError("residue_subset must be 'specific' or 'nucleotide'")

    df = df.loc[df["condition"] == condition].copy()

    ranked = (
        df.groupby("residue_label", as_index=False)[metric]
        .mean()
        .sort_values(metric, ascending=False)
        .head(top_n)
    )

    output_path = output_dir / f"top_{top_n}_{metric}_{suffix}_{condition}.png"

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(ranked["residue_label"], ranked[metric])
    ax.set_title(f"Top {top_n} {suffix} residues by {metric.replace('_', ' ')} ({condition})")
    ax.set_xlabel("Residue")
    ax.set_ylabel(metric.replace("_", " ").title())
    ax.set_ylim(bottom=0)
    ax.tick_params(axis="x", rotation=45)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()

    return output_path