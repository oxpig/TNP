"""The flag thresholds for TNP's metrics.

These are the single definition of TNP's thresholds, which were calibrated on the
clinical-stage nanobodies described in the TNP paper.  They are used to flag
each metric, to draw the threshold lines on the web application's plots, and to
tabulate the thresholds on the web application's home page.
"""

import math
from dataclasses import dataclass

#: TNP's metrics, by the abbreviations used in its output.
METRICS = {
    "L": "Total CDR Length",
    "L3": "CDR3 Length",
    "C": "CDR3 Compactness",
    "PSH": "Patches of Surface Hydrophobicity",
    "PPC": "Patches of Positive Charge",
    "PNC": "Patches of Negative Charge",
}


@dataclass(frozen=True)
class Thresholds:
    """The flag regions of a metric.

    As their names suggest, `red_below` and `red_above` are strict, and
    `amber_up_to` and `amber_from` are inclusive.  A value is flagged red if it
    is below `red_below` or above `red_above`; amber if it is up to and
    including `amber_up_to`, or from `amber_from` up; and green otherwise.  So
    each amber region is closed, and a value on any threshold is flagged amber.
    This keeps the extremes of the calibration set, which `red_below` and
    `red_above` are, out of the red regions.

    A metric with no lower flag regions has `red_below` and `amber_up_to` set to
    -∞.
    """

    red_below: int | float
    amber_up_to: int | float
    amber_from: int | float
    red_above: int | float

    @property
    def two_sided(self):
        """Whether the metric is flagged for low values as well as high."""
        return self.red_below > -math.inf

    def flag(self, value):
        """Return the flag colour, "red", "amber" or "green", for a value.

        A metric that could not be calculated (e.g. compactness, if anchor
        residues are missing) is flagged "grey", rather than falling through every
        comparison below to green.
        """
        if value is None or math.isnan(value):
            return "grey"
        if value < self.red_below or value > self.red_above:
            return "red"
        if value <= self.amber_up_to or value >= self.amber_from:
            return "amber"
        return "green"


# These unrounded values are those of the plot templates' threshold lines, and of
# the web application's table of thresholds.  Until they were defined here, flags
# were assigned with rounded values instead:
#     - L 37,
#     - L3 21,
#     - C 0.56, 0.81, 1.57 & 1.61,
#     - PSH 73.40, 79.59, 126.83 & 155.47,
#     - PPC 0.39 & 1.18,
#     - PNC 1.47 & 1.88.
THRESHOLDS = {
    # CDR lengths are whole numbers of residues, so their thresholds are too.  The
    # calibrated amber thresholds for long loops, 37.25 and 21.25, are equivalent
    # to 38 and 22.
    "L": Thresholds(red_below=20, amber_up_to=24, amber_from=38, red_above=39),
    "L3": Thresholds(red_below=5, amber_up_to=8, amber_from=22, red_above=23),
    "C": Thresholds(
        red_below=0.5634161539313097,
        amber_up_to=0.8145579042527696,
        amber_from=1.5691348616692684,
        red_above=1.613227748114984,
    ),
    "PSH": Thresholds(
        red_below=73.4039,
        amber_up_to=79.594375,
        amber_from=126.82575,
        red_above=155.4699,
    ),
    "PPC": Thresholds(
        red_below=-math.inf,
        amber_up_to=-math.inf,
        amber_from=0.394125,
        red_above=1.1853,
    ),
    "PNC": Thresholds(
        red_below=-math.inf,
        amber_up_to=-math.inf,
        amber_from=1.470925,
        red_above=1.8838,
    ),
}


def assign_flag(metric, value):
    """Return the flag colour for a metric's value: "red", "amber" or "green", or
    "grey" if it could not be calculated.
    """
    return THRESHOLDS[metric].flag(value)
