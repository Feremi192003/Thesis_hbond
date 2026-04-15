from __future__ import annotations

import argparse
import time
from pathlib import Path

from .pipeline import run_pipeline

""" Min run point, this is what parses user arguments, translates into pipeline options, does the analysis and prints/graphs outputs.
"""

def build_parser() -> argparse.ArgumentParser:
    """Create the command-line argument parser for the package.

    Returns:
        argparse.ArgumentParser: Parser configured with I/O, optional flags to contorl if you want plotting, residue metrics or parallelism
    """
    parser = argparse.ArgumentParser(
        description="Hydrogen-bond interaction analysis pipeline for p53 isoform CSV files."
    )

    #This is needed bc it checks for input directory containing raw input files post CPPTRAJ
    parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing *_DNA_RES_INTERACT.csv files",
    )
    #Output directory where stuff should be written.
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("results"),
        help="Directory where summary CSVs and plots will be written",
    )
    #Number of workers being used if we are using parallel mode.
    parser.add_argument(
        "-w",
        "--workers",
        type=int,
        default=4,
        help="Number of worker processes for parallel analysis",
    )
    #Allow for process based exe. across multiple files.
    parser.add_argument(
        "--parallel",
        action="store_true",
        help="Enable process-based parallel execution",
    )
    #skip plotting.
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip plot generation",
    )
    #skip row by row annotation tables.
    parser.add_argument(
        "--skip-row-annotations",
        action="store_true",
        help="Do not generate row-by-row annotation output",
    )
    #Skip residue metrics and only give file summaries.
    parser.add_argument(
        "--skip-residue-metrics",
        action="store_true",
        help="Do not generate residue-level metric outputs",
    )
    #For user convenience, skips most parser.add_arguments, mostly used for benchmarking case study in this thesis.
    parser.add_argument(
        "--metrics-only",
        action="store_true",
        help="Fast benchmarking mode: disables plots, row annotations, and residue metrics",
    )

    return parser


def main() -> None:
    """Parse arguments, run pipeline, print short summary.
    """
    parser = build_parser()
    args = parser.parse_args()

    #Translate flags into boolean.
    make_plots = not args.no_plots
    include_row_annotations = not args.skip_row_annotations
    include_residue_metrics = not args.skip_residue_metrics

    #If user only wants benchmarking disable everything that could potentially slow down main-file-level options.
    if args.metrics_only:
        make_plots = False
        include_row_annotations = False
        include_residue_metrics = False

    #Overall timer for performance engineering
    start = time.perf_counter()

    outputs = run_pipeline(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        make_plots=make_plots,
        workers=args.workers,
        use_parallel=args.parallel,
        include_row_annotations=include_row_annotations,
        include_residue_metrics=include_residue_metrics,
        write_profile=True,
    )

    end = time.perf_counter()

    print("\nPipeline completed. Generated files:")
    for label, path in outputs.items():
        print(f"- {label}: {path}")

    print(f"\nTotal runtime (main.py wrapper): {end - start:.4f} seconds")


if __name__ == "__main__":
    main()