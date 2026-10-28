from theraprofnano.structure.process_pdb import get_modelling_details


def test_modelling_details(clinical_model):
    model = clinical_model("seq10")
    details = get_modelling_details(str(model), type="nanobody", save=True)

    # DSSP ran, and annotated every residue.
    assert details["dssp_message"] == ""
    assert len(details["dssp_properties"]["H"]) == 127
    assert (model.parent / "modelling_details.jsonp").is_file()
    assert details["imgt_CDRH_ranges"]
    assert details["n_possible_seq_liabilities"] >= 0
