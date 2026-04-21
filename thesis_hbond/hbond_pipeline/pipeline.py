from __future__ import annotations

"""Top-down runs hbond workflow, this module puts together every other module into a huge run.
"""

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import math
import os
import time

import pandas as pd

from .aggregate import (
    aggregate_metrics,
    aggregate_residue_metrics,
    build_condition_tables,
    condition_comparison_table,
)
from .config import CANONICAL_TARGET_RESIDUES, DEFAULT_PATTERN
from .metrics import compute_file_metrics
from .parsing import list_csv_files, parse_filename_metadata, read_interaction_csv
from .plots import (
    plot_metric_by_isoform,
    plot_residue_counts,
    plot_residue_metric_by_isoform,
    plot_residue_metric_heatmap,
    plot_top_residue_metric,
)


def _process_single_csv(
    csv_file: Path,
    include_row_annotations: bool,
    include_residue_metrics: bool,
) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    """Read one csv file, compute metrics, return resulting tables.
    This helper exists so the same logic can be used for serial and parallel modes.
    Returns:
        tuple[dict, pd.DataFrame, pd.DataFrame]:
    """
    metadata = parse_filename_metadata(csv_file)
    df = read_interaction_csv(csv_file)

    metrics_row, annotated_df, residue_df = compute_file_metrics(
        df=df,
        metadata=metadata,
        canonical_target_residues=CANONICAL_TARGET_RESIDUES,
        include_row_annotations=include_row_annotations,
        include_residue_metrics=include_residue_metrics,
    )

    return metrics_row, annotated_df, residue_df


def _infer_chunksize(n_files: int, workers: int) -> int:
    """Choose executor chunk size for parallel processing.
    """
    if workers <= 1:
        return 1
    return max(1, math.ceil(n_files / (workers * 4)))


def run_pipeline(
    input_dir: str | Path,
    output_dir: str | Path,
    pattern: str = DEFAULT_PATTERN,
    make_plots: bool = True,
    workers: int = 1,
    use_parallel: bool = False,
    include_row_annotations: bool = True,
    include_residue_metrics: bool = True,
    write_profile: bool = True,
) -> dict[str, Path]:
    """Run complete pipeline

    Args:
    input_dir, output_dir:
        Input location for raw interaction CSV files and output location for all
        generated tables/plots.
    pattern:
        Recursive used to discover CSV files.
    make_plots:
        Whether to generate plots in addition to CSV tables.
    workers:
        Maximum number of worker processes to use when ``use_parallel`` is true.
    use_parallel:
        Whether to process files with ``ProcessPoolExecutor``.
    include_row_annotations:
        Whether to generate per-row annotation tables.
    include_residue_metrics:
        Whether to generate residue-level tables and plots.
    write_profile:
        Whether to export a CSV timing breakdown for the pipeline.

    Returns:
        dict[str, Path]: mapping from output names to the files written on disk.
    """
    total_start = time.perf_counter()

    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    #Store per stage runtimes for performance engineering and analysis.
    timings: dict[str, float] = {}

    # File discovery
    t0 = time.perf_counter()
    csv_files = list_csv_files(input_path, pattern=pattern)
    t1 = time.perf_counter()
    timings["file_discovery_s"] = t1 - t0

    if not csv_files:
        raise FileNotFoundError(f"No CSV files matching '{pattern}' found in {input_path}")

    #Cap worker count on maxiple number of CPU cores acailable.
    max_cpu = os.cpu_count() or 1
    workers = max(1, min(int(workers), max_cpu))

    mode = "parallel" if use_parallel else "serial"
    print(f"Found CSV files: {len(csv_files)}")
    print(f"Execution mode: {mode}")
    print(f"Workers used: {workers}")
    print(f"Row annotations enabled: {include_row_annotations}")
    print(f"Residue metrics enabled: {include_residue_metrics}")
    print(f"Plot generation enabled: {make_plots}")
    print(f"File discovery: {timings['file_discovery_s']:.4f} s")

    metrics_rows: list[dict] = []
    annotated_row_tables: list[pd.DataFrame] = []
    residue_metric_tables: list[pd.DataFrame] = []

    # Metric computation
    t2 = time.perf_counter()

    if not use_parallel:
        for csv_file in csv_files:
            metrics_row, annotated_df, residue_df = _process_single_csv(
                csv_file=csv_file,
                include_row_annotations=include_row_annotations,
                include_residue_metrics=include_residue_metrics,
            )
            metrics_rows.append(metrics_row)
            if include_row_annotations:
                annotated_row_tables.append(annotated_df)
            if include_residue_metrics:
                residue_metric_tables.append(residue_df)
    else:
        chunksize = _infer_chunksize(len(csv_files), workers)
        print(f"Parallel chunksize: {chunksize}")

        with ProcessPoolExecutor(max_workers=workers) as executor:
            results = executor.map(
                _process_single_csv,
                csv_files,
                [include_row_annotations] * len(csv_files),
                [include_residue_metrics] * len(csv_files),
                chunksize=chunksize,
            )
            for metrics_row, annotated_df, residue_df in results:
                metrics_rows.append(metrics_row)
                if include_row_annotations:
                    annotated_row_tables.append(annotated_df)
                if include_residue_metrics:
                    residue_metric_tables.append(residue_df)

    t3 = time.perf_counter()
    timings["metric_computation_s"] = t3 - t2
    print(f"Metric computation: {timings['metric_computation_s']:.4f} s")

    # 3. Aggregation
    t4 = time.perf_counter()

    file_metrics_df = pd.DataFrame(metrics_rows).sort_values(
        ["isoform", "condition", "replicate", "filename"]
    )

    aggregated_df = aggregate_metrics(file_metrics_df)
    comparison_df = condition_comparison_table(aggregated_df)
    wt_table_df, y220c_table_df = build_condition_tables(aggregated_df)

    if include_row_annotations and annotated_row_tables:
        row_annotations_df = pd.concat(annotated_row_tables, ignore_index=True)
    else:
        row_annotations_df = pd.DataFrame()

    if include_residue_metrics and residue_metric_tables:
        file_residue_metrics_df = pd.concat(residue_metric_tables, ignore_index=True)
        aggregated_residue_df = aggregate_residue_metrics(file_residue_metrics_df)

        #these tables are used for targeted exports and plots
        specific_residue_df = aggregated_residue_df.loc[
            aggregated_residue_df["is_specific_residue"]
        ].copy()

        nucleotide_residue_df = aggregated_residue_df.loc[
            aggregated_residue_df["is_nucleotide_residue"]
        ].copy()
    else:
        file_residue_metrics_df = pd.DataFrame()
        aggregated_residue_df = pd.DataFrame()
        specific_residue_df = pd.DataFrame()
        nucleotide_residue_df = pd.DataFrame()

    t5 = time.perf_counter()
    timings["aggregation_s"] = t5 - t4
    print(f"Aggregation/tables: {timings['aggregation_s']:.4f} s")

    # Output file paths 
    file_metrics_csv = output_path / "file_level_metrics.csv"
    aggregated_csv = output_path / "isoform_condition_averages.csv"
    comparison_csv = output_path / "wt_y220c_comparison.csv"
    wt_table_csv = output_path / "wt_isoform_metrics.csv"
    y220c_table_csv = output_path / "y220c_isoform_metrics.csv"
    summary_totals_csv = output_path / "summary_totals.csv"

    row_annotations_csv = output_path / "row_annotations.csv"
    file_residue_metrics_csv = output_path / "file_level_residue_metrics.csv"
    aggregated_residue_csv = output_path / "residue_condition_isoform_averages.csv"
    specific_residue_csv = output_path / "specific_residue_averages.csv"
    nucleotide_residue_csv = output_path / "nucleotide_residue_averages.csv"
    timings_csv = output_path / "pipeline_timings.csv"

    # Write outputs
    t6 = time.perf_counter()

    file_metrics_df.to_csv(file_metrics_csv, index=False)
    aggregated_df.to_csv(aggregated_csv, index=False)
    comparison_df.to_csv(comparison_csv, index=False)
    wt_table_df.to_csv(wt_table_csv, index=False)
    y220c_table_df.to_csv(y220c_table_csv, index=False)
    comparison_df.to_csv(summary_totals_csv, index=False)

    if include_row_annotations:
        row_annotations_df.to_csv(row_annotations_csv, index=False)

    if include_residue_metrics:
        file_residue_metrics_df.to_csv(file_residue_metrics_csv, index=False)
        aggregated_residue_df.to_csv(aggregated_residue_csv, index=False)
        specific_residue_df.to_csv(specific_residue_csv, index=False)
        nucleotide_residue_df.to_csv(nucleotide_residue_csv, index=False)

    t7 = time.perf_counter()
    timings["write_outputs_s"] = t7 - t6
    print(f"Writing outputs: {timings['write_outputs_s']:.4f} s")

    outputs: dict[str, Path] = {
        "file_level_metrics": file_metrics_csv,
        "isoform_condition_averages": aggregated_csv,
        "wt_y220c_comparison": comparison_csv,
        "wt_isoform_metrics": wt_table_csv,
        "y220c_isoform_metrics": y220c_table_csv,
        "summary_totals": summary_totals_csv,
    }

    if include_row_annotations:
        outputs["row_annotations"] = row_annotations_csv

    if include_residue_metrics:
        outputs["file_level_residue_metrics"] = file_residue_metrics_csv
        outputs["residue_condition_isoform_averages"] = aggregated_residue_csv
        outputs["specific_residue_averages"] = specific_residue_csv
        outputs["nucleotide_residue_averages"] = nucleotide_residue_csv

    # Plot generation if requested.
    t8 = time.perf_counter()

    if make_plots:
        outputs["sis_plot"] = plot_metric_by_isoform(aggregated_df, file_metrics_df, "sis", output_path)
        outputs["nucleotide_sis_plot"] = plot_metric_by_isoform(
            aggregated_df, file_metrics_df, "nucleotide_sis", output_path
        )
        outputs["ehbs_plot"] = plot_metric_by_isoform(aggregated_df, file_metrics_df, "ehbs", output_path)
        outputs["trt_plot"] = plot_metric_by_isoform(aggregated_df, file_metrics_df, "trt", output_path)

        residue_plots = plot_residue_counts(aggregated_df, file_metrics_df, output_path)
        for path in residue_plots:
            outputs[path.stem] = path

        if include_residue_metrics and not aggregated_residue_df.empty:
            outputs["specific_sis_heatmap"] = plot_residue_metric_heatmap(
                aggregated_residue_df,
                metric="sis_contribution",
                output_dir=output_path,
                residue_subset="specific",
            )
            outputs["nucleotide_sis_heatmap"] = plot_residue_metric_heatmap(
                aggregated_residue_df,
                metric="nucleotide_sis_contribution",
                output_dir=output_path,
                residue_subset="nucleotide",
            )
            outputs["nucleotide_occupancy_heatmap"] = plot_residue_metric_heatmap(
                aggregated_residue_df,
                metric="total_occupancy",
                output_dir=output_path,
                residue_subset="nucleotide",
            )
            outputs["nucleotide_contacts_heatmap"] = plot_residue_metric_heatmap(
                aggregated_residue_df,
                metric="unique_dna_contacts",
                output_dir=output_path,
                residue_subset="nucleotide",
            )

            outputs["nucleotide_occupancy_plot"] = plot_residue_metric_by_isoform(
                aggregated_residue_df,
                metric="total_occupancy",
                output_dir=output_path,
                residue_subset="nucleotide",
            )
            outputs["nucleotide_contacts_plot"] = plot_residue_metric_by_isoform(
                aggregated_residue_df,
                metric="unique_dna_contacts",
                output_dir=output_path,
                residue_subset="nucleotide",
            )
            outputs["nucleotide_sis_contribution_plot"] = plot_residue_metric_by_isoform(
                aggregated_residue_df,
                metric="nucleotide_sis_contribution",
                output_dir=output_path,
                residue_subset="nucleotide",
            )

            outputs["top_specific_sis_wt"] = plot_top_residue_metric(
                aggregated_residue_df,
                metric="sis_contribution",
                output_dir=output_path,
                condition="WT",
                residue_subset="specific",
                top_n=8,
            )
            outputs["top_specific_sis_y220c"] = plot_top_residue_metric(
                aggregated_residue_df,
                metric="sis_contribution",
                output_dir=output_path,
                condition="Y220C",
                residue_subset="specific",
                top_n=8,
            )

    t9 = time.perf_counter()
    timings["plot_generation_s"] = t9 - t8
    print(f"Plot generation: {timings['plot_generation_s']:.4f} s")

    #Record total runtime and export a timing profile.
    total_end = time.perf_counter()
    timings["total_runtime_s"] = total_end - total_start
    print(f"Total runtime: {timings['total_runtime_s']:.4f} s")

    if write_profile:
        timings_df = pd.DataFrame(
            [{"stage": key, "seconds": value} for key, value in timings.items()]
        )
        timings_df.to_csv(timings_csv, index=False)
        outputs["pipeline_timings"] = timings_csv

    return outputs