from theraprofnano.thresholds import THRESHOLDS


def set_threshold_lines(graph, metric):
    """Draw a metric's threshold lines on a plot template at its thresholds.

    The templates' first data series holds the x positions of the threshold
    lines: red, amber, amber and red for two-sided metrics, in columns 0, 2, 3
    and 4; and amber and red for one-sided metrics, in columns 0 and 2.
    """
    thresholds = THRESHOLDS[metric]
    if thresholds.two_sided:
        lines = {
            0: thresholds.red_below,
            2: thresholds.amber_up_to,
            3: thresholds.amber_from,
            4: thresholds.red_above,
        }
    else:
        lines = {0: thresholds.amber_from, 2: thresholds.red_above}
    for row in graph["data"]["data01"]:
        for column, x in lines.items():
            row[column] = x
