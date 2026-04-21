from __future__ import annotations
""" This is the module that standardizes how files are discovered, 
    extracts metadata from their names, 
    converts canonical res numbers into iso-specific and reads raw CSVS into pandas.
"""
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd

from .config import COLUMN_NAMES, FILENAME_PATTERNS, ISOFORM_OFFSETS

# Residue labels are expected to look like "ARG_248".
RESIDUE_RE = re.compile(r"^([A-Z]{3})_(\d+)$")


@dataclass(frozen=True)
class FileMetadata:
    """Container for metadata parsed from an input filename, in order:
    filepath:
        Full path to the input CSV file.
    filename:
        Basename of the file including extension.
    stem:
        Filename without extension.
    condition:
        Experimental condition, normalized to WT or Y220C.
    isoform:
        Isoform number extracted from the filename.
    replicate:
        Replicate number if present, otherwise None.
    offset:
        Sequence truncation offset used to map canonical residues onto the
        isoform's residue numbering.

    """
    filepath: Path
    filename: str
    stem: str
    condition: str
    isoform: int
    replicate: int | None
    offset: int


def parse_filename_metadata(filepath: str | Path) -> FileMetadata:
    """Extract condition, iso, rep, offset info from filename.
    """
    path = Path(filepath)
    name = path.name
    stem = path.stem

    #match independently, so if parser fails it can send an error message saying what is missing from file name.
    condition_match = re.search(FILENAME_PATTERNS["condition"], stem)
    iso_match = re.search(FILENAME_PATTERNS["isoform"], stem)
    rep_match = re.search(FILENAME_PATTERNS["replicate"], stem)

    if not condition_match:
        raise ValueError(f"Could not parse condition from filename: {name}")
    if not iso_match:
        raise ValueError(f"Could not parse isoform from filename: {name}")

    raw_condition = condition_match.group(1)
    condition = "WT" if raw_condition == "WT" else "Y220C"
    isoform = int(iso_match.group(1))
    replicate = int(rep_match.group(1)) if rep_match else None
    offset = ISOFORM_OFFSETS.get(isoform, 0)

    return FileMetadata(
        filepath=path,
        filename=name,
        stem=stem,
        condition=condition,
        isoform=isoform,
        replicate=replicate,
        offset=offset,
    )


def adjust_residue_for_isoform(canonical_residue: str, offset: int) -> str | None:
    """Converts canonical residues into the numbering of isoform, 
    returns None when substracting the offset if it would be result in 0 or a negative index.

    """
    match = RESIDUE_RE.match(canonical_residue)
    if not match:
        raise ValueError(f"Invalid residue label: {canonical_residue}")

    aa, residue_number_text = match.groups()
    adjusted_number = int(residue_number_text) - offset
    if adjusted_number <= 0:
        return None

    return f"{aa}_{adjusted_number}"


def adjusted_residue_map(canonical_residues: Iterable[str], offset: int) -> dict[str, str]:
    """Maps canonical residues to iso-adjusted labels for one offset. Residues that are not in that offset are ommited from dictionary.
    """
    mapping: dict[str, str] = {}
    for residue in canonical_residues:
        adjusted = adjust_residue_for_isoform(residue, offset)
        if adjusted is not None:
            mapping[residue] = adjusted
    return mapping


def read_interaction_csv(filepath: str | Path) -> pd.DataFrame:
    """read one CSV into a pandas DtaFrame.
    
    The function inspects the first non-comment data line to decide if the file is comma-spearated or whitespace.
    
    Reads the data using the columns defined in config.py

    Returns:
        pd.DataFrame.
    """
    path = Path(filepath)

    first_data_line = None
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                first_data_line = stripped
                break

    if first_data_line is None:
        raise ValueError(f"No readable data lines found in file: {path}")

    #Choose the delimiter dynamically, some files may be comma separated and others whitespaced.
    sep = "," if "," in first_data_line else r"\s+"

    df = pd.read_csv(
        path,
        sep=sep,
        engine="python",
        header=None,
        names=COLUMN_NAMES,
        comment="#",
        dtype=str,
    )

    #Make sure all columns exist even if input file is shorter than expected.
    for col in COLUMN_NAMES:
        if col not in df.columns:
            df[col] = ""

    df = df.fillna("")
    return df


def list_csv_files(input_path: str | Path, pattern: str = "*_DNA_RES_INTERACT.csv") -> list[Path]:
    """
    Discover matching interaction csv files under input_path recursively.
    """
    root = Path(input_path)
    return sorted(root.rglob(pattern))