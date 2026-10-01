import json
import subprocess

import pytest

SEQUENCE = (
    "QVKLQESGAELARPGASVKLSCKASGYTFTNYWMQWVKQRPGQGLDWIGAIYPGDGNTRYTHKFKGKATLTADKSSST"
    "AYMQLSSLASEDSGVYYCARGEGNYAWFAYWGQGTTVTVSS"
)


@pytest.mark.slow
def test_single_sequence(tmp_path):
    """Run TNP on a sequence, including modelling it with NanoBodyBuilder2."""
    subprocess.run(
        ["TNP", "--seq", SEQUENCE, "--output", str(tmp_path), "--web"], check=True
    )
    results = json.loads((tmp_path / "TNP_Results_SingleSeqEntry_Nb1.json").read_text())
    nb1 = results["Nb1"]
    assert nb1["Total CDR Length"] == 28
    assert nb1["CDR3 Length"] == 12
    assert set(nb1["Flags"]) == {"L", "L3", "C", "PSH", "PPC", "PNC"}
    assert (tmp_path / "PSH_D3_Template.json").is_file()
    assert (tmp_path / "Final_Models" / "Nb1_NanoBodyBuilder2_Model_Details.json").is_file()
