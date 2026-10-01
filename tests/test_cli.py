import subprocess

import numpy as np
import pytest

from theraprofnano.cli import assign_flag


def test_help():
    result = subprocess.run(["TNP", "--help"], capture_output=True, text=True, check=True)
    assert "--seq" in result.stdout


@pytest.mark.parametrize(
    "metric, value, flag",
    [
        ("L", 28, "green"),
        ("L", 22, "amber"),
        ("L", 40, "red"),
        ("C", 0.7, "amber"),
        ("PSH", 75.0, "amber"),
        ("PSH", 79.9, "green"),
        ("PPC", 0.5, "amber"),
        ("PNC", 2.0, "red"),
    ],
)
def test_assign_flag(metric, value, flag):
    assert assign_flag(metric, value) == flag


@pytest.mark.parametrize("metric", ["L", "L3", "C", "PSH", "PPC", "PNC"])
def test_uncomputable_metric_is_flagged_red(metric):
    assert assign_flag(metric, np.nan) == "red"
