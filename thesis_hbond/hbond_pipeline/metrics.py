from __future__ import annotations

"""This is the core metric computation for a single CSV.
    It contains the single-pass streaming algorithm that converts raw rows into file-level metrics, row annotations and residue-level contribution tables.
"""

from collections import defaultdict
from dataclasses import asdict
from typing import Any
import re

import pandas as pd

from .config import NUCLEOTIDE_TARGET_RESIDUES
from .parsing import FileMetadata, adjusted_residue_map

#Regex to extract lables like ARG_248 from cells like ARG_248@NH1
PROTEIN_RESIDUE_RE = re.compile(r"\b([A-Z]{3}_\d+)@")


#Mapping for rows returned by df.itertuples, avoids repeated string-based column lookups inside the main loop.
# itertuples(index=False, name=None) positions for config.COLUMN_NAMES:
# 0 dna_atom
# 1 protein_atom_1
# 2 protein_atom_2
# 3 frames
# 4 frac
# 5 avgdist
# 6 avgang
DNA_IDX = 0
PROT1_IDX = 1
PROT2_IDX = 2
FRAC_IDX = 4


def _safe_float(value: Any) -> float:
    """Converts values to floats, return 0.0 if fails. Function is mainly used to handle frac, which may contian something other than a number.
    """
   
    try:
        out = pd.to_numeric(value, errors="coerce")
        if pd.isna(out):
            return 0.0
        return float(out)
    except Exception:
        return 0.0


def _extract_protein_residues_from_strings(protein_atom_1: str, protein_atom_2: str) -> list[str]:
    """Extract unique residue labels from the protein atom columns.

    Takes raw string fields from an hbond interaction row, 

    Returns sorted unique res lables found in two protein atom strings, we are essentially avoiding double-counting of residues.
    """
    residues: set[str] = set()

    #search protein fields because either column may contain a residue that means something to the interaction.
    if protein_atom_1:
        residues.update(PROTEIN_RESIDUE_RE.findall(protein_atom_1))
    if protein_atom_2:
        residues.update(PROTEIN_RESIDUE_RE.findall(protein_atom_2))

    return sorted(residues)


def compute_file_metrics(
    df: pd.DataFrame,
    metadata: FileMetadata,
    canonical_target_residues: list[str],
    include_row_annotations: bool = True,
    include_residue_metrics: bool = True,
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
    """
    Single-pass metric computation for one interaction CSV.

    Returns:
        result: file-level metric dictionary
        annotated_df: optional row-level annotation table
        residue_df: optional residue-level contribution table
    """
    #Maps canonical residue lables onto the numbering scheme for the speciifc isoform its working on by substracting to the appropiate truncation offset.
    residue_mapping = adjusted_residue_map(canonical_target_residues, metadata.offset)
    nucleotide_mapping = adjusted_residue_map(NUCLEOTIDE_TARGET_RESIDUES, metadata.offset)

    #Reverse lookups so adjusted residue strings in the CSV can be mapped back to canonical labels.
    adjusted_to_canonical = {
        adjusted: canonical for canonical, adjusted in residue_mapping.items()
    }
    adjusted_to_nucleotide = {
        adjusted: canonical for canonical, adjusted in nucleotide_mapping.items()
    }

    adjusted_specific_residues = tuple(residue_mapping.values())
    adjusted_nucleotide_residues = tuple(nucleotide_mapping.values())
    
    #Track row-unique counts for canonical target residues.
    residue_counts = {canonical: 0 for canonical in residue_mapping}

    #Residue-level accumulators used only when residue metrics are enabled.
    occupancy_by_residue: defaultdict[str, float] = defaultdict(float)
    dna_contacts_by_residue: defaultdict[str, set[str]] = defaultdict(set)
    sis_contribution_by_residue: defaultdict[str, int] = defaultdict(int)
    nucleotide_sis_contribution_by_residue: defaultdict[str, int] = defaultdict(int)

    annotated_rows: list[dict[str, Any]] = []

    #Initialize file-level metrics.
    sis = 0
    nucleotide_sis = 0
    ehbs = len(df)

    frac_values = pd.to_numeric(df["frac"], errors="coerce").fillna(0.0)
    trt = float(frac_values.sum())
    
    
    #Iterate ONCE through the table so all outputs can be built from the same pass over data.
    for row_number, row in enumerate(df.itertuples(index=False, name=None), start=1):
        dna_atom = row[DNA_IDX] if row[DNA_IDX] else ""
        protein_atom_1 = row[PROT1_IDX] if row[PROT1_IDX] else ""
        protein_atom_2 = row[PROT2_IDX] if row[PROT2_IDX] else ""
        frac_value = _safe_float(row[FRAC_IDX])

        #Combine relevant strings into one searchable record for checks below.
        combined_text = f"{dna_atom} {protein_atom_1} {protein_atom_2}"

        matched_adjusted = [res for res in adjusted_specific_residues if res in combined_text]
        matched_nucleotide_adjusted = [res for res in adjusted_nucleotide_residues if res in combined_text]

        matched_canonical = [adjusted_to_canonical[res] for res in matched_adjusted]
        matched_nucleotide_canonical = [adjusted_to_nucleotide[res] for res in matched_nucleotide_adjusted]

        row_has_match = len(matched_adjusted) > 0
        row_has_nucleotide_match = len(matched_nucleotide_adjusted) > 0


        #SIS Counts if a row has a target residue, residue_counts tracks which canonical targets contributed to that file-level value.
        if row_has_match:
            sis += 1
            for adjusted_residue in matched_adjusted:
                canonical_residue = adjusted_to_canonical[adjusted_residue]
                residue_counts[canonical_residue] += 1

        #Same idea for nuc-sis
        if row_has_nucleotide_match:
            nucleotide_sis += 1
        #Parse out residues appearing in this row so their TRT can be acculumated.
        if include_residue_metrics:
            protein_residues = _extract_protein_residues_from_strings(protein_atom_1, protein_atom_2)

            for residue_label in protein_residues:
                occupancy_by_residue[residue_label] += frac_value

                if dna_atom:
                    dna_contacts_by_residue[residue_label].add(dna_atom)

                if residue_label in adjusted_to_canonical:
                    sis_contribution_by_residue[residue_label] += 1

                if residue_label in adjusted_to_nucleotide:
                    nucleotide_sis_contribution_by_residue[residue_label] += 1

        if include_row_annotations:
            #Save explanation of what matched in this row so user can check after.
            annotated_rows.append(
                {
                    "filename": metadata.filename,
                    "condition": metadata.condition,
                    "isoform": metadata.isoform,
                    "replicate": metadata.replicate,
                    "row_number": row_number,
                    "row_has_target_match": row_has_match,
                    "row_has_nucleotide_match": row_has_nucleotide_match,
                    "matched_canonical_residues": ";".join(matched_canonical),
                    "matched_adjusted_residues": ";".join(matched_adjusted),
                    "matched_nucleotide_canonical_residues": ";".join(matched_nucleotide_canonical),
                    "matched_nucleotide_adjusted_residues": ";".join(matched_nucleotide_adjusted),
                }
            )

    #Put the file level metrics and metadata into one row.
    result: dict[str, Any] = {
        **asdict(metadata),
        "sis": sis,
        "nucleotide_sis": nucleotide_sis,
        "ehbs": ehbs,
        "trt": trt,
    }

    #add one count column per residue so the layer can build residue specific bar charts automatically.
    for canonical_residue in canonical_target_residues:
        result[f"count_{canonical_residue}"] = residue_counts.get(canonical_residue, 0)

    #Include every residue athat appeared in residue-level accumulators so the output table is complete.
    if include_residue_metrics:
        all_residues = sorted(
            set(occupancy_by_residue)
            | set(dna_contacts_by_residue)
            | set(sis_contribution_by_residue)
            | set(nucleotide_sis_contribution_by_residue)
        )

        residue_rows: list[dict[str, Any]] = []
        for residue_label in all_residues:
            canonical_specific = adjusted_to_canonical.get(residue_label, "")
            canonical_nucleotide = adjusted_to_nucleotide.get(residue_label, "")
            #Prefer canonical lables for tracked residues, else keep original residue lables.
            canonical_label = canonical_specific or canonical_nucleotide or residue_label 


            residue_rows.append(
                {
                    "filename": metadata.filename,
                    "condition": metadata.condition,
                    "isoform": metadata.isoform,
                    "replicate": metadata.replicate,
                    "residue_label": canonical_label,
                    "canonical_specific_residue": canonical_specific,
                    "canonical_nucleotide_residue": canonical_nucleotide,
                    "is_specific_residue": bool(canonical_specific),
                    "is_nucleotide_residue": bool(canonical_nucleotide),
                    "total_occupancy": occupancy_by_residue.get(residue_label, 0.0),
                    "unique_dna_contacts": len(dna_contacts_by_residue.get(residue_label, set())),
                    "sis_contribution": sis_contribution_by_residue.get(residue_label, 0),
                    "nucleotide_sis_contribution": nucleotide_sis_contribution_by_residue.get(residue_label, 0),
                }
            )

        residue_df = pd.DataFrame(residue_rows)
    else:
        #Return an empty table with expected colums even if disabled, so code runs properly regardless.
        residue_df = pd.DataFrame(
            columns=[
                "filename",
                "condition",
                "isoform",
                "replicate",
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

    if include_row_annotations:
        annotated_df = pd.DataFrame(annotated_rows)
    else:
        #Same as comment above, so callers of this function can still check or inspect.
        annotated_df = pd.DataFrame(
            columns=[
                "filename",
                "condition",
                "isoform",
                "replicate",
                "row_number",
                "row_has_target_match",
                "row_has_nucleotide_match",
                "matched_canonical_residues",
                "matched_adjusted_residues",
                "matched_nucleotide_canonical_residues",
                "matched_nucleotide_adjusted_residues",
            ]
        )

    return result, annotated_df, residue_df