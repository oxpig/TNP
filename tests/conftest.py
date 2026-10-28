import shutil
from pathlib import Path

import pytest

CLINICAL_MODELS = (
    Path(__file__).resolve().parent.parent
    / "paper"
    / "paper_data"
    / "vhh_clinical_set_models"
)


@pytest.fixture
def clinical_model(tmp_path, monkeypatch):
    """Copy a clinical-stage nanobody model into a temporary working directory.

    TNP writes its outputs beside its inputs and into the working directory.
    """
    monkeypatch.chdir(tmp_path)

    def copy(name):
        return Path(shutil.copy(CLINICAL_MODELS / f"{name}_NanoBodyBuilder2_Model.pdb", tmp_path))

    return copy


def pytest_addoption(parser):
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)
