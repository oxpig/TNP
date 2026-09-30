"""Residue side-chain accessibility, computed with FreeSASA.

This reproduces the output of the `psa` program (version 2.0, from JOY) that TNP
used to call, on which TNP's metric thresholds were calibrated.  `psa` reports
*contact* surface areas: the area of each atom's van der Waals sphere that a
1.4 Å probe can touch.  FreeSASA computes solvent-accessible areas, from which
the contact area of each atom follows exactly, by scaling by (r / (r + probe))².

With NACCESS atomic radii and the Lee & Richards algorithm, this matches `psa`'s
side-chain areas to within its printed precision (0.01 Å²) on the 36 clinical-
stage nanobody models used to calibrate TNP.
"""

import freesasa
from Bio.PDB import PDBParser

PROBE_RADIUS = 1.4

# NACCESS atomic radii (Å), as tabulated by FreeSASA.
_C_ALIPHATIC, _C_AROMATIC, _N, _O, _S = 1.87, 1.76, 1.65, 1.40, 1.85
_BACKBONE_RADII = {"N": _N, "CA": _C_ALIPHATIC, "C": _C_AROMATIC, "O": _O,
                   "OXT": _O, "CB": _C_ALIPHATIC}
_SIDE_CHAIN_RADII = {
    "ARG": {"CG": _C_ALIPHATIC, "CD": _C_ALIPHATIC, "NE": _N, "CZ": _C_AROMATIC,
            "NH1": _N, "NH2": _N},
    "ASN": {"CG": _C_AROMATIC, "OD1": _O, "ND2": _N},
    "ASP": {"CG": _C_AROMATIC, "OD1": _O, "OD2": _O},
    "CYS": {"SG": _S},
    "GLN": {"CG": _C_ALIPHATIC, "CD": _C_AROMATIC, "OE1": _O, "NE2": _N},
    "GLU": {"CG": _C_ALIPHATIC, "CD": _C_AROMATIC, "OE1": _O, "OE2": _O},
    "HIS": {"CG": _C_AROMATIC, "ND1": _N, "CD2": _C_AROMATIC, "CE1": _C_AROMATIC,
            "NE2": _N},
    "ILE": {"CG1": _C_ALIPHATIC, "CG2": _C_ALIPHATIC, "CD1": _C_ALIPHATIC},
    "LEU": {"CG": _C_ALIPHATIC, "CD1": _C_ALIPHATIC, "CD2": _C_ALIPHATIC},
    "LYS": {"CG": _C_ALIPHATIC, "CD": _C_ALIPHATIC, "CE": _C_ALIPHATIC,
            "NZ": 1.50},
    "MET": {"CG": _C_ALIPHATIC, "SD": _S, "CE": _C_ALIPHATIC},
    "PHE": {"CG": _C_AROMATIC, "CD1": _C_AROMATIC, "CD2": _C_AROMATIC,
            "CE1": _C_AROMATIC, "CE2": _C_AROMATIC, "CZ": _C_AROMATIC},
    "PRO": {"CG": _C_ALIPHATIC, "CD": _C_ALIPHATIC},
    "SER": {"OG": _O},
    "THR": {"OG1": _O, "CG2": _C_ALIPHATIC},
    "TRP": {"CG": _C_AROMATIC, "CD1": _C_AROMATIC, "CD2": _C_AROMATIC,
            "NE1": _N, "CE2": _C_AROMATIC, "CE3": _C_AROMATIC,
            "CZ2": _C_AROMATIC, "CZ3": _C_AROMATIC, "CH2": _C_AROMATIC},
    "TYR": {"CG": _C_AROMATIC, "CD1": _C_AROMATIC, "CD2": _C_AROMATIC,
            "CE1": _C_AROMATIC, "CE2": _C_AROMATIC, "CZ": _C_AROMATIC, "OH": _O},
    "VAL": {"CG1": _C_ALIPHATIC, "CG2": _C_ALIPHATIC},
}
# Fallbacks, by element, for atoms of non-standard residues and hetero groups.
_ELEMENT_RADII = {"C": 1.80, "N": 1.60, "O": 1.40, "S": 1.85}

# `psa` counts every atom but the backbone N, C and O as side chain, so that
# CA, and the C-terminal OXT, contribute to a residue's side-chain area.
_BACKBONE_ATOMS = {"N", "C", "O"}

# `psa`'s side-chain contact areas (Å²) for each residue type in an extended
# Ala-X-Ala tripeptide, against which its relative accessibilities are
# expressed.  These are recovered from `psa`'s own output, as sum / per cent.
REFERENCE_SIDE_CHAIN_AREAS = {
    "ALA": 23.160, "ARG": 62.660, "ASN": 30.932, "ASP": 29.399, "CYS": 31.950,
    "GLN": 41.672, "GLU": 38.610, "GLY": 10.840, "HIS": 45.612, "ILE": 45.711,
    "LEU": 46.608, "LYS": 51.529, "MET": 51.944, "PHE": 51.536, "PRO": 39.281,
    "SER": 23.620, "THR": 31.930, "TRP": 65.080, "TYR": 52.912, "VAL": 37.887,
}


def _radius(residue_name, atom_name, element):
    """Return the NACCESS radius of an atom."""
    radius = _SIDE_CHAIN_RADII.get(residue_name, {}).get(atom_name)
    if radius is None and residue_name in REFERENCE_SIDE_CHAIN_AREAS:
        radius = _BACKBONE_RADII.get(atom_name)
    if radius is None:
        radius = _ELEMENT_RADII.get(element, 1.80)
    return radius


def side_chain_accessibility(pdb_file):
    """Calculate the side-chain contact surface area of each residue.

    Hydrogen atoms and waters are ignored, as by `psa`.

    Returns two dictionaries, keyed by chain identifier and then by residue
    identifier (the residue number followed by any insertion code, e.g. "112A"):
    the relative side-chain accessibility of each residue, in per cent, and its
    absolute side-chain contact area, in Å².  These are rounded to the precision
    that `psa` printed, so that TNP's thresholds behave as they did with `psa`.
    Residues of unknown type are given only an absolute area.
    """
    structure = PDBParser(QUIET=True).get_structure("model", pdb_file)
    model = next(structure.get_models())

    coordinates, radii, atoms = [], [], []
    for atom in model.get_atoms():
        residue = atom.get_parent()
        if atom.element in {"H", "D"} or residue.get_resname() == "HOH":
            continue
        radius = _radius(residue.get_resname(), atom.get_id(), atom.element)
        coordinates.extend(float(x) for x in atom.coord)
        radii.append(radius)
        atoms.append(atom)

    parameters = freesasa.Parameters(
        {"algorithm": freesasa.LeeRichards, "probe-radius": PROBE_RADIUS}
    )
    result = freesasa.calcCoord(coordinates, radii, parameters)

    areas = {}
    for i, atom in enumerate(atoms):
        if atom.get_id() in _BACKBONE_ATOMS:
            continue
        residue = atom.get_parent()
        chain = residue.get_parent().id
        _, number, insertion_code = residue.id
        key = (chain, f"{number}{insertion_code.strip()}", residue.get_resname())
        contact_area = result.atomArea(i) * (radii[i] / (radii[i] + PROBE_RADIUS)) ** 2
        areas[key] = areas.get(key, 0.0) + contact_area

    relative, absolute = {}, {}
    for (chain, residue_id, residue_name), area in areas.items():
        absolute.setdefault(chain, {})[residue_id] = round(area, 2)
        if residue_name in REFERENCE_SIDE_CHAIN_AREAS:
            reference = REFERENCE_SIDE_CHAIN_AREAS[residue_name]
            relative.setdefault(chain, {})[residue_id] = round(
                100 * area / reference, 1
            )
    return relative, absolute
