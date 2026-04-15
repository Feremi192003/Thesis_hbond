from __future__ import annotations
"""configurating the residues that must be looked for, CHANGE IF YOU'RE NOT USING P53!!!!

Stores canonical residue targets, isoform offsets, filename parsing.
"""

#Canonical p53 residues between p53-DNA contacts with DNA=1TUP
CANONICAL_TARGET_RESIDUES = [
    "SER_241",
    "ARG_283",
    "CYS_277",
    "ALA_276",
    "ARG_280",
    "ARG_273",
    "LYS_120",
    "ARG_248",
]

#Those residues that specifically bind to the DNA Amino Acid and not the DNA Backbone
NUCLEOTIDE_TARGET_RESIDUES = [
    "LYS_120",
    "CYS_277",
    "ARG_280",
]

#Isoform groups truncate different parts of the N-Terminal so offseting here for that.
ISOFORM_OFFSETS = {
    1: 0,
    2: 0,
    3: 0,
    4: 39,
    5: 39,
    6: 39,
    7: 132,
    8: 132,
    9: 132,
    10: 159,
    11: 159,
    12: 159,
}


#RegEx  to extract condition, isoform and replicate values frm filenames
FILENAME_PATTERNS = {
    "condition": r"^(WT|MUT)",
    "isoform": r"(?:^|[_\s.-])i(\d{1,2})(?:\b|[_\s.-])",
    "replicate": r"(?:^|[_\s.-])r(\d{1,2})(?:\b|[_\s.-])",
}

#What column order should look like from CPPTRAJ Files.
COLUMN_NAMES = [
    "dna_atom",
    "protein_atom_1",
    "protein_atom_2",
    "frames",
    "frac",
    "avgdist",
    "avgang",
]
#Recursiing to find CSV files.
DEFAULT_PATTERN = "*_DNA_RES_INTERACT.csv"