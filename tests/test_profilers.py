"""Regression tests on clinical-stage nanobody models from the TNP paper."""

import pytest
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

from theraprofnano.CDR_Profiler.CDR3_Conf_Assigner import main_compactness
from theraprofnano.CDR_Profiler.CDR_Assigner import main as cdr_profile
from theraprofnano.Hydrophobicity_and_Charge_Profiler.Hydrophobicity_and_Charge_Assigner import (  # noqa: E501
    CreateAnnotation,
)

# Calculated with `psa`, before it was replaced with FreeSASA.
PATCHES = {
    "seq1": (79.3304, 0.0515, 0.0),
    "seq10": (88.7075, 0.1082, 1.4503),
    "seq18": (79.0911, 1.1123, 0.0925),
}
CDR_LENGTHS = {
    "seq1": {"cdrh1": 8, "cdrh2": 7, "cdrh3": 5},
    "seq10": {"cdrh1": 10, "cdrh2": 8, "cdrh3": 18},
    "seq18": {"cdrh1": 8, "cdrh2": 7, "cdrh3": 15},
}
CDR3_REACH = {"seq1": 8.8761, "seq10": 12.1866, "seq18": 14.8088}


@pytest.mark.parametrize("name", PATCHES)
def test_surface_patches(clinical_model, name):
    stats = CreateAnnotation(0, 7.4, str(clinical_model(name)), "IG", "imgt", verbose=False)
    psh, ppc, pnc = PATCHES[name]
    assert stats[0]["Patch_Hydrophob_CDR"] == pytest.approx(psh, abs=1e-4)
    assert stats[0]["Patch_Pos_Charge_CDR"] == pytest.approx(ppc, abs=1e-4)
    assert stats[0]["Patch_Neg_Charge_CDR"] == pytest.approx(pnc, abs=1e-4)


@pytest.mark.parametrize("name", CDR_LENGTHS)
def test_cdr_lengths(clinical_model, name):
    model = clinical_model(name)
    residues = PDBParser(QUIET=True).get_structure(name, model).get_residues()
    sequence = "".join(seq1(r.get_resname()) for r in residues)
    _, lengths = cdr_profile(f"{name}_H", sequence, "H", str(model.parent), verbose=False)
    assert lengths["imgt"] == CDR_LENGTHS[name]


@pytest.mark.parametrize("name", CDR3_REACH)
def test_cdr3_reach(clinical_model, name):
    reach = main_compactness(str(clinical_model(name)), "imgt", verbose=False)
    assert reach == pytest.approx(CDR3_REACH[name], abs=1e-4)
