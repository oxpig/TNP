"""
A module to deal with region annotations for different schemes. To replace the regions tuples.
"""

regions_in_chothia = {
    "L1": {
        "kabat": (24, 34),
        "chothia": (24, 34),
        "contact": (30, 36),
        "north": (24, 34),
    },
    "L2": {
        "kabat": (50, 56),
        "chothia": (50, 56),
        "contact": (46, 55),
        "north": (49, 56),
    },
    "L3": {
        "kabat": (89, 97),
        "chothia": (89, 97),
        "contact": (89, 96),
        "north": (89, 97),
    },
    "H1": {
        "kabat": (31, 35),
        "chothia": (26, 32),
        "contact": (30, 35),
        "north": (23, 35),
    },
    "H2": {
        "kabat": (50, 65),
        "chothia": (52, 56),
        "contact": (47, 58),
        "north": (50, 58),
    },
    "H3": {
        "kabat": (95, 102),
        "chothia": (95, 102),
        "contact": (93, 101),
        "north": (93, 102),
    },
}


_regions = {"imgt": {}}
_regions["imgt"]["L"] = _regions["imgt"]["H"] = (
    "11111111111111111111111111222222222222333333333333333334444444444555555555555555555555555555555555555555666666666666677777777777"
)
_regions["kabat"] = {}
_regions["kabat"]["L"] = (
    "11111111111111111111111222222222222222223333333333333334444444444444455555555555555555555555555555555555666666666666677777777777"
)
_regions["kabat"]["H"] = (
    "11111111111111111111111111111111111222223333333333333344444444444444444444555555555555555555555555555555556666666666677777777777"
)
_regions["chothia"] = {}
_regions["chothia"]["L"] = (
    "11111111111111111111111222222222222222223333333333333334444444444444455555555555555555555555555555555555666666666666677777777777"
)
_regions["chothia"]["H"] = (
    "11111111111111111111111111222222222223333333333333333333444444445555555555555555555555555555555555555555556666666666677777777777"
)
_regions["contact"] = {}
_regions["contact"]["L"] = (
    "11111111111111111111111111111111111222222233333333344444444444444444555555555555555555555555555555555555666666666666777777777777"
)
_regions["contact"]["H"] = (
    "11111111111111111111111111111122222222223333333333344444444444444455555555555555555555555555555555555555666666666666777777777777"
)
_regions["north"] = {}
_regions["north"]["L"] = (
    "11111111111111111111111222222222222222223333333333333344444444444444455555555555555555555555555555555555666666666666677777777777"
)
_regions["north"]["H"] = (
    "11111111111111111111111222222222222222223333333333333344444444444455555555555555555555555555555555555555666666666666677777777777"
)


# For internal use only. These are not direct conversions and are handled heuristically.
_index_to_imgt_state = {
    ("chothia", "H"): {
        1: 0,
        2: 1,
        3: 2,
        4: 3,
        5: 4,
        6: 6,
        7: 7,
        8: 8,
        9: 9,
        10: 10,
        11: 11,
        12: 12,
        13: 13,
        14: 14,
        15: 15,
        16: 16,
        17: 17,
        18: 18,
        19: 19,
        20: 20,
        21: 21,
        22: 22,
        23: 23,
        24: 24,
        25: 25,
        26: 26,
        27: 27,
        28: 28,
        29: 29,
        30: 30,
        31: 35,
        32: 36,
        33: 37,
        34: 38,
        35: 39,
        36: 40,
        37: 41,
        38: 42,
        39: 43,
        40: 44,
        41: 45,
        42: 46,
        43: 47,
        44: 48,
        45: 49,
        46: 50,
        47: 51,
        48: 52,
        49: 53,
        50: 54,
        51: 55,
        52: 59,
        53: 60,
        54: 61,
        55: 62,
        56: 63,
        57: 64,
        58: 65,
        59: 66,
        60: 67,
        61: 68,
        62: 69,
        63: 70,
        64: 72,
        65: 73,
        66: 74,
        67: 75,
        68: 76,
        69: 77,
        70: 78,
        71: 79,
        72: 80,
        73: 81,
        74: 82,
        75: 83,
        76: 84,
        77: 85,
        78: 86,
        79: 87,
        80: 88,
        81: 89,
        82: 93,
        83: 94,
        84: 95,
        85: 96,
        86: 97,
        87: 98,
        88: 99,
        89: 100,
        90: 101,
        91: 102,
        92: 103,
        93: 104,
        94: 105,
        95: 106,
        96: 107,
        97: 108,
        98: 109,
        99: 110,
        100: 114,
        101: 115,
        102: 116,
        103: 117,
        104: 118,
        105: 119,
        106: 120,
        107: 121,
        108: 122,
        109: 123,
        110: 124,
        111: 125,
        112: 126,
        113: 127,
    },
    ("kabat", "H"): {
        1: 0,
        2: 1,
        3: 2,
        4: 3,
        5: 4,
        6: 6,
        7: 7,
        8: 8,
        9: 9,
        10: 10,
        11: 11,
        12: 12,
        13: 13,
        14: 14,
        15: 15,
        16: 16,
        17: 17,
        18: 18,
        19: 19,
        20: 20,
        21: 21,
        22: 22,
        23: 23,
        24: 24,
        25: 25,
        26: 26,
        27: 27,
        28: 28,
        29: 29,
        30: 30,
        31: 31,
        32: 32,
        33: 33,
        34: 34,
        35: 35,
        36: 40,
        37: 41,
        38: 42,
        39: 43,
        40: 44,
        41: 45,
        42: 46,
        43: 47,
        44: 48,
        45: 49,
        46: 50,
        47: 51,
        48: 52,
        49: 53,
        50: 54,
        51: 55,
        52: 59,
        53: 60,
        54: 61,
        55: 62,
        56: 63,
        57: 64,
        58: 65,
        59: 66,
        60: 67,
        61: 68,
        62: 69,
        63: 70,
        64: 72,
        65: 73,
        66: 74,
        67: 75,
        68: 76,
        69: 77,
        70: 78,
        71: 79,
        72: 80,
        73: 81,
        74: 82,
        75: 83,
        76: 84,
        77: 85,
        78: 86,
        79: 87,
        80: 88,
        81: 89,
        82: 93,
        83: 94,
        84: 95,
        85: 96,
        86: 97,
        87: 98,
        88: 99,
        89: 100,
        90: 101,
        91: 102,
        92: 103,
        93: 104,
        94: 105,
        95: 106,
        96: 107,
        97: 108,
        98: 109,
        99: 110,
        100: 114,
        101: 115,
        102: 116,
        103: 117,
        104: 118,
        105: 119,
        106: 120,
        107: 121,
        108: 122,
        109: 123,
        110: 124,
        111: 125,
        112: 126,
        113: 127,
    },
    ("imgt", "H"): {
        1: 0,
        2: 1,
        3: 2,
        4: 3,
        5: 4,
        6: 5,
        7: 6,
        8: 7,
        9: 8,
        10: 9,
        11: 10,
        12: 11,
        13: 12,
        14: 13,
        15: 14,
        16: 15,
        17: 16,
        18: 17,
        19: 18,
        20: 19,
        21: 20,
        22: 21,
        23: 22,
        24: 23,
        25: 24,
        26: 25,
        27: 26,
        28: 27,
        29: 28,
        30: 29,
        31: 30,
        32: 31,
        33: 32,
        34: 33,
        35: 34,
        36: 35,
        37: 36,
        38: 37,
        39: 38,
        40: 39,
        41: 40,
        42: 41,
        43: 42,
        44: 43,
        45: 44,
        46: 45,
        47: 46,
        48: 47,
        49: 48,
        50: 49,
        51: 50,
        52: 51,
        53: 52,
        54: 53,
        55: 54,
        56: 55,
        57: 56,
        58: 57,
        59: 58,
        60: 59,
        61: 60,
        62: 61,
        63: 62,
        64: 63,
        65: 64,
        66: 65,
        67: 66,
        68: 67,
        69: 68,
        70: 69,
        71: 70,
        72: 71,
        73: 72,
        74: 73,
        75: 74,
        76: 75,
        77: 76,
        78: 77,
        79: 78,
        80: 79,
        81: 80,
        82: 81,
        83: 82,
        84: 83,
        85: 84,
        86: 85,
        87: 86,
        88: 87,
        89: 88,
        90: 89,
        91: 90,
        92: 91,
        93: 92,
        94: 93,
        95: 94,
        96: 95,
        97: 96,
        98: 97,
        99: 98,
        100: 99,
        101: 100,
        102: 101,
        103: 102,
        104: 103,
        105: 104,
        106: 105,
        107: 106,
        108: 107,
        109: 108,
        110: 109,
        111: 110,
        112: 111,
        113: 112,
        114: 113,
        115: 114,
        116: 115,
        117: 116,
        118: 117,
        119: 118,
        120: 119,
        121: 120,
        122: 121,
        123: 122,
        124: 123,
        125: 124,
        126: 125,
        127: 126,
        128: 127,
    },
    ("chothia", "L"): {
        1: 0,
        2: 1,
        3: 2,
        4: 3,
        5: 4,
        6: 5,
        7: 6,
        8: 7,
        9: 8,
        10: 9,
        11: 10,
        12: 11,
        13: 12,
        14: 13,
        15: 14,
        16: 15,
        17: 16,
        18: 17,
        19: 18,
        20: 19,
        21: 20,
        22: 21,
        23: 22,
        24: 23,
        25: 24,
        26: 25,
        27: 26,
        28: 27,
        29: 28,
        30: 35,
        31: 36,
        32: 37,
        33: 38,
        34: 39,
        35: 40,
        36: 41,
        37: 42,
        38: 43,
        39: 44,
        40: 45,
        41: 46,
        42: 47,
        43: 48,
        44: 49,
        45: 50,
        46: 51,
        47: 52,
        48: 53,
        49: 54,
        50: 55,
        51: 56,
        52: 57,
        53: 65,
        54: 66,
        55: 67,
        56: 68,
        57: 69,
        58: 70,
        59: 72,
        60: 73,
        61: 74,
        62: 75,
        63: 76,
        64: 77,
        65: 78,
        66: 81,
        67: 82,
        68: 83,
        69: 84,
        70: 85,
        71: 86,
        72: 87,
        73: 88,
        74: 89,
        75: 90,
        76: 91,
        77: 92,
        78: 93,
        79: 94,
        80: 95,
        81: 96,
        82: 97,
        83: 98,
        84: 99,
        85: 100,
        86: 101,
        87: 102,
        88: 103,
        89: 104,
        90: 105,
        91: 106,
        92: 107,
        93: 108,
        94: 109,
        95: 114,
        96: 115,
        97: 116,
        98: 117,
        99: 118,
        100: 119,
        101: 120,
        102: 121,
        103: 122,
        104: 123,
        105: 124,
        106: 125,
        107: 126,
        108: 127,
    },
    ("martin", "H"): {
        1: 0,
        2: 1,
        3: 2,
        4: 3,
        5: 4,
        6: 5,
        7: 6,
        8: 8,
        9: 9,
        10: 10,
        11: 11,
        12: 12,
        13: 13,
        14: 14,
        15: 15,
        16: 16,
        17: 17,
        18: 18,
        19: 19,
        20: 20,
        21: 21,
        22: 22,
        23: 23,
        24: 24,
        25: 25,
        26: 26,
        27: 27,
        28: 28,
        29: 29,
        30: 30,
        31: 35,
        32: 36,
        33: 37,
        34: 38,
        35: 39,
        36: 40,
        37: 41,
        38: 42,
        39: 43,
        40: 44,
        41: 45,
        42: 46,
        43: 47,
        44: 48,
        45: 49,
        46: 50,
        47: 51,
        48: 52,
        49: 53,
        50: 54,
        51: 55,
        52: 59,
        53: 60,
        54: 61,
        55: 62,
        56: 63,
        57: 64,
        58: 65,
        59: 66,
        60: 67,
        61: 68,
        62: 69,
        63: 70,
        64: 72,
        65: 73,
        66: 74,
        67: 75,
        68: 76,
        69: 77,
        70: 78,
        71: 79,
        72: 83,
        73: 84,
        74: 85,
        75: 86,
        76: 87,
        77: 88,
        78: 89,
        79: 90,
        80: 91,
        81: 92,
        82: 93,
        83: 94,
        84: 95,
        85: 96,
        86: 97,
        87: 98,
        88: 99,
        89: 100,
        90: 101,
        91: 102,
        92: 103,
        93: 104,
        94: 105,
        95: 106,
        96: 107,
        97: 108,
        98: 109,
        99: 110,
        100: 114,
        101: 115,
        102: 116,
        103: 117,
        104: 118,
        105: 119,
        106: 120,
        107: 121,
        108: 122,
        109: 123,
        110: 124,
        111: 125,
        112: 126,
        113: 127,
    },
    ("kabat", "L"): {
        1: 0,
        2: 1,
        3: 2,
        4: 3,
        5: 4,
        6: 5,
        7: 6,
        8: 7,
        9: 8,
        10: 9,
        11: 10,
        12: 11,
        13: 12,
        14: 13,
        15: 14,
        16: 15,
        17: 16,
        18: 17,
        19: 18,
        20: 19,
        21: 20,
        22: 21,
        23: 22,
        24: 23,
        25: 24,
        26: 25,
        27: 32,
        28: 33,
        29: 34,
        30: 35,
        31: 36,
        32: 37,
        33: 38,
        34: 39,
        35: 40,
        36: 41,
        37: 42,
        38: 43,
        39: 44,
        40: 45,
        41: 46,
        42: 47,
        43: 48,
        44: 49,
        45: 50,
        46: 51,
        47: 52,
        48: 53,
        49: 54,
        50: 55,
        51: 56,
        52: 57,
        53: 65,
        54: 66,
        55: 67,
        56: 68,
        57: 69,
        58: 70,
        59: 72,
        60: 73,
        61: 74,
        62: 75,
        63: 76,
        64: 77,
        65: 78,
        66: 81,
        67: 82,
        68: 83,
        69: 84,
        70: 85,
        71: 86,
        72: 87,
        73: 88,
        74: 89,
        75: 90,
        76: 91,
        77: 92,
        78: 93,
        79: 94,
        80: 95,
        81: 96,
        82: 97,
        83: 98,
        84: 99,
        85: 100,
        86: 101,
        87: 102,
        88: 103,
        89: 104,
        90: 105,
        91: 106,
        92: 107,
        93: 108,
        94: 109,
        95: 114,
        96: 115,
        97: 116,
        98: 117,
        99: 118,
        100: 119,
        101: 120,
        102: 121,
        103: 122,
        104: 123,
        105: 124,
        106: 125,
        107: 126,
        108: 127,
    },
    ("imgt", "L"): {
        1: 0,
        2: 1,
        3: 2,
        4: 3,
        5: 4,
        6: 5,
        7: 6,
        8: 7,
        9: 8,
        10: 9,
        11: 10,
        12: 11,
        13: 12,
        14: 13,
        15: 14,
        16: 15,
        17: 16,
        18: 17,
        19: 18,
        20: 19,
        21: 20,
        22: 21,
        23: 22,
        24: 23,
        25: 24,
        26: 25,
        27: 26,
        28: 27,
        29: 28,
        30: 29,
        31: 30,
        32: 31,
        33: 32,
        34: 33,
        35: 34,
        36: 35,
        37: 36,
        38: 37,
        39: 38,
        40: 39,
        41: 40,
        42: 41,
        43: 42,
        44: 43,
        45: 44,
        46: 45,
        47: 46,
        48: 47,
        49: 48,
        50: 49,
        51: 50,
        52: 51,
        53: 52,
        54: 53,
        55: 54,
        56: 55,
        57: 56,
        58: 57,
        59: 58,
        60: 59,
        61: 60,
        62: 61,
        63: 62,
        64: 63,
        65: 64,
        66: 65,
        67: 66,
        68: 67,
        69: 68,
        70: 69,
        71: 70,
        72: 71,
        73: 72,
        74: 73,
        75: 74,
        76: 75,
        77: 76,
        78: 77,
        79: 78,
        80: 79,
        81: 80,
        82: 81,
        83: 82,
        84: 83,
        85: 84,
        86: 85,
        87: 86,
        88: 87,
        89: 88,
        90: 89,
        91: 90,
        92: 91,
        93: 92,
        94: 93,
        95: 94,
        96: 95,
        97: 96,
        98: 97,
        99: 98,
        100: 99,
        101: 100,
        102: 101,
        103: 102,
        104: 103,
        105: 104,
        106: 105,
        107: 106,
        108: 107,
        109: 108,
        110: 109,
        111: 110,
        112: 111,
        113: 112,
        114: 113,
        115: 114,
        116: 115,
        117: 116,
        118: 117,
        119: 118,
        120: 119,
        121: 120,
        122: 121,
        123: 122,
        124: 123,
        125: 124,
        126: 125,
        127: 126,
        128: 127,
    },
    ("martin", "L"): {
        1: 0,
        2: 1,
        3: 2,
        4: 3,
        5: 4,
        6: 5,
        7: 6,
        8: 7,
        9: 8,
        10: 9,
        11: 10,
        12: 11,
        13: 12,
        14: 13,
        15: 14,
        16: 15,
        17: 16,
        18: 17,
        19: 18,
        20: 19,
        21: 20,
        22: 21,
        23: 22,
        24: 23,
        25: 24,
        26: 25,
        27: 26,
        28: 27,
        29: 28,
        30: 35,
        31: 36,
        32: 37,
        33: 38,
        34: 39,
        35: 40,
        36: 41,
        37: 42,
        38: 43,
        39: 44,
        40: 45,
        41: 46,
        42: 47,
        43: 48,
        44: 49,
        45: 50,
        46: 51,
        47: 52,
        48: 53,
        49: 54,
        50: 55,
        51: 56,
        52: 57,
        53: 65,
        54: 66,
        55: 67,
        56: 68,
        57: 69,
        58: 70,
        59: 72,
        60: 73,
        61: 74,
        62: 75,
        63: 76,
        64: 77,
        65: 78,
        66: 81,
        67: 82,
        68: 83,
        69: 84,
        70: 85,
        71: 86,
        72: 87,
        73: 88,
        74: 89,
        75: 90,
        76: 91,
        77: 92,
        78: 93,
        79: 94,
        80: 95,
        81: 96,
        82: 97,
        83: 98,
        84: 99,
        85: 100,
        86: 101,
        87: 102,
        88: 103,
        89: 104,
        90: 105,
        91: 106,
        92: 107,
        93: 108,
        94: 109,
        95: 114,
        96: 115,
        97: 116,
        98: 117,
        99: 118,
        100: 119,
        101: 120,
        102: 121,
        103: 122,
        104: 123,
        105: 124,
        106: 125,
        107: 126,
        108: 127,
    },
}

wolfguy_indexdiv50_to_region = {
    "H": ["fwh1", "cdrh1", "fwh2", "cdrh2", "fwh3", "cdrh3", "fwh4"],
    "L": ["fwl1", "cdrl1", "fwl2", "cdrl2", "fwl3", "cdrl3", "fwl4"],
}

_reg_one2three = {
    "1": "fw%s1",
    "2": "cdr%s1",
    "3": "fw%s2",
    "4": "cdr%s2",
    "5": "fw%s3",
    "6": "cdr%s3",
    "7": "fw%s4",
}


def get_region(position, chain, numbering_scheme="chothia", definition="chothia"):
    """
    Get the region in which the position belongs given the chain, numbering scheme and definition.

    **Note** this function does not know about insertions on the sequence. Therefore, it will get the region annotation
    wrong when using non-equivalent scheme-definitions.

    To get around this please use the annotate_regions function which implements heuristics to get the definition correct
    in the scheme.

    """
    index, insertion = position
    chain = chain.upper()

    # Horrible exception cases revolving around the kabat scheme/definition and cdr h1
    if definition == "kabat":
        if (
            numbering_scheme == "kabat" and chain == "H" and 31 <= index < 36
        ):  # Kabat scheme kabat definition.
            if index == 35:
                if insertion in " AB":  # Position 31 to 35B
                    return "cdrh1"
                else:
                    return "fwh2"  # 31C would be framework.
            else:
                return "cdrh1"
    if numbering_scheme == "kabat":  # Kabat numbering, chothia or imgt definitions.
        if definition == "chothia" and chain == "H" and 33 <= index < 36:
            return "fwh2"
        elif definition == "imgt" and chain == "H" and 34 <= index < 36:
            return "fwh2"

    if numbering_scheme == "wolfguy" or definition == "wolfguy":
        assert definition == "wolfguy" and numbering_scheme == "wolfguy", (
            "The wolfguy numbering scheme must be used with the wolfguy CDR definition"
        )
        if chain == "H":
            if index > 411:
                return ""
            r = int(index / 50) - 2
        elif chain == "L":
            if index > 810:
                return ""
            r = int(index / 50) - 10
        try:
            return wolfguy_indexdiv50_to_region[chain][r]
        except IndexError:
            return ""

    try:
        return (
            _reg_one2three[
                _regions[definition][chain][
                    _index_to_imgt_state[(numbering_scheme, chain)][index]
                ]
            ]
            % chain.lower()
        )
    except KeyError:
        return "?"


def annotate_regions(
    numbered_sequence, chain, numbering_scheme="chothia", definition="chothia"
):
    """
    Given a numbered sequence annotate which region each residue belongs to.

    The numbering scheme can be one chothia, kabat, imgt or martin
    The definition can be chothia, kabat, imgt, north or contact.

    Contact definition cannot be used with the kabat numbering scheme.

    If possible, use the corresponding numbering scheme and definition.

    This function automates the heuristics recognise different definitions in each scheme. However,
    some of the conversions are non-trivial.

    """
    # In some cases there is not a direct equivalence between numbering schemes. Therefore, a CDR definition
    # may not be easily translatable from scheme to scheme. This is the case for:
    #   - H1 Kabat   definition in the IMGT    scheme   # 1
    #   - H1 Chothia definition in the Kabat   scheme   # 2
    #   - H1 IMGT    definition in the Kabat   scheme   # 3
    #   - H1 Chothia definition in the Kabat   scheme   # 4
    #   - L2 IMGT    definition in the Kabat   scheme   # 4
    #   - L2 IMGT    definition in the Chothia scheme   # 5
    #   - L2 IMGT    definition in the Martin  scheme   # 6

    # Please stop using the Kabat scheme and the Kabat definition, it is a pain the arse.

    # Below are heuristics to allow the conversion. They have been developed so that the CDR sequence extracted
    # in all schemes is the same as the CDR sequence in the scheme the definition was originally defined. e.g.
    # Chothia for Chothia, Kabat for Kabat, IMGT for IMGT. Note that the Contact definition *cannot* be defined
    # in Kabat numbering.

    # Find which additional positions should be considered in the CDR, given the sequence.
    additional_positions = {}
    excluded_positions = {}
    c = chain.lower()

    numdict = dict(numbered_sequence)

    cdr_acceptors = {
        1: Accept(numbering_scheme=numbering_scheme, definition=definition),
        2: Accept(numbering_scheme=numbering_scheme, definition=definition),
        3: Accept(numbering_scheme=numbering_scheme, definition=definition),
    }

    cdr_acceptors[1].set_regions(["cdr%s1" % c])
    cdr_acceptors[2].set_regions(["cdr%s2" % c])
    cdr_acceptors[3].set_regions(["cdr%s3" % c])

    # Heavy chain
    if chain == "H":
        if (
            numbering_scheme == "imgt" and definition == "kabat"
        ):  # IMGT scheme / Kabat definition
            # 1
            # Count the positions 31-34 inclusive. These become insertions on 35 in Kabat.
            ins = 0
            for i in range(31, 35):
                if (i, " ") in numdict:
                    if numdict[(i, " ")] != "-":
                        ins += 1
            for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
                if (33, a) in numdict:
                    if numdict[(33, a)] != "-":
                        ins += 1
                else:
                    break
            # If insertions would occur on 35 in the Kabat scheme, extend the IMGT numbering back until you hit 31. Then put on the insertions at 32
            if ins:  # Add postions backwards based on the number of insertions starting with H35.
                # IMGT insertions go on 33
                cdr_acceptors[1].add_positions(
                    (
                        [(35, " "), (34, " "), (33, " "), (32, " ")]
                        + [(33, a) for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"]
                    )[:ins],
                    "H",
                )
                # If there more than two insertions we have to start removing positions from the definition. See # 4 below for more of an explanation
                if ins > 2:
                    cdr_acceptors[1].exclude_positions(
                        (
                            [
                                (40, " "),
                                (39, " "),
                                (38, " "),
                                (37, " "),
                                (36, " "),
                                (35, " "),
                                (34, " "),
                            ]
                            + [(33, a) for a in "FGHIJKLMNOPQRSTUVWXYZ"]
                        )[: ins - 2],
                        "H",
                    )
        elif numbering_scheme == "kabat" and definition in [
            "chothia",
            "imgt",
        ]:  # Kabat scheme / Chothia or IMGT definition
            # 2, 3
            # Count the insertions on 35. These happen outside the CDR definitions of IMGT or Chothia.
            # They therefore are missed by a straight conversion. The equivalence is dependent on the
            # number of insertions at this position.

            # Count the number of insertions on H35.
            ins = 0
            for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
                if (35, a) in numdict:
                    if numdict[(35, a)] != "-":
                        ins += 1
                else:
                    break
            if ins:  # Add postions based on the number of insertions starting with H33 for the Chothia definition. H34 for the IMGT definition.
                if definition == "chothia":
                    cdr_acceptors[1].add_positions(
                        (
                            [(33, " "), (34, " ")]
                            + [(35, a) for a in " ABCDEFGHIJKLMNOPQRSTUVWXYZ"]
                        )[:ins],
                        "H",
                    )
                elif definition == "imgt":
                    cdr_acceptors[1].add_positions(
                        (
                            [(34, " ")]
                            + [(35, a) for a in " ABCDEFGHIJKLMNOPQRSTUVWXYZ"]
                        )[:ins],
                        "H",
                    )
        elif numbering_scheme in ["chothia", "martin"] and definition == "kabat":
            # 4
            # Count the insertions on 31. If there are more than two, exclude back from 35.
            # e.g. if we have the below. Kabat and Chothia numbered Kabat CDRs
            # [((31, ' '), 'P'), ((32, ' '), 'A'), ((33, ' '), 'P'), ((34, ' '), 'E'), ((35, ' '), 'H'), ((35, 'A'), 'F'), ((35, 'B'), 'I')]
            # [((31, ' '), 'P'), ((31, 'A'), 'A'), ((31, 'B'), 'P'), ((31, 'C'), 'E'), ((32, ' '), 'H'), ((33, ' '), 'F'), ((34, ' '), 'I')]
            # Kabat definition will have a max of length 7 as insertions are on 35 and the definition ends at 35B.
            # In Chothia numbering the insertions are at 31. Hence when more than 2 insertions occur Chothia 35 falls off the end. Chothia 34
            # becomes the end of the Kabat CDR. Similar thing happens with IMGT.
            # Count the number of insertions on H31.
            ins = 0
            for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
                if (31, a) in numdict:
                    if numdict[(31, a)] != "-":
                        ins += 1
                else:
                    break
            if ins > 2:
                cdr_acceptors[1].exclude_positions(
                    (
                        [(35, " "), (34, " "), (33, " "), (32, " ")]
                        + [(31, a) for a in "GHIJKLMNOPQRSTUVWXYZ"]
                    )[: ins - 2],
                    "H",
                )

    # Light chain
    if chain == "L":
        # 5,6,7
        if (
            numbering_scheme in ["kabat", "chothia", "martin"] and definition == "imgt"
        ):  # Kabat-like schemes and IMGT definition
            # Count the number of insertions on L54.
            ins = 0
            for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
                if (54, a) in numdict:
                    if numdict[(54, a)] != "-":
                        ins += 1
                else:
                    break
            if (
                ins
            ):  # Add positions based on the number of insertions starting with L53.
                # IMGT definition with no insertions ends and 52 in Kabat/Chothia/Martin.
                extensions = [(53, " ")] + [
                    (54, a) for a in " ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                ]
                cdr_acceptors[2].add_positions(
                    ([(53, " ")] + [(54, a) for a in " ABCDEFGHIJKLMNOPQRSTUVWXYZ"])[
                        :ins
                    ],
                    "L",
                )

    fw_regions = ["fw%s1" % c, "fw%s2" % c, "fw%s3" % c, "fw%s4" % c]
    fw_region = "fw%s1" % c
    region_annotations = []
    cterm = max(_index_to_imgt_state[(numbering_scheme, chain)].keys())
    for r, a in numbered_sequence:
        if cdr_acceptors[1].accept(r, chain):
            region_annotations.append((r, a, "cdr%s1" % c))
            fw_region = "fw%s2" % c
        elif cdr_acceptors[2].accept(r, chain):
            region_annotations.append((r, a, "cdr%s2" % c))
            fw_region = "fw%s3" % c
        elif cdr_acceptors[3].accept(r, chain):
            region_annotations.append((r, a, "cdr%s3" % c))
            fw_region = "fw%s4" % c
        elif (
            r[0] <= cterm
        ):  # Anything out of the variable region is not assigned a region i.e. ''
            region_annotations.append((r, a, fw_region))
        else:
            region_annotations.append((r, a, ""))

    return region_annotations


class Accept:
    """
    A class to select which positions should be compared.
    """

    _defined_regions = [
        "fwh1",
        "fwh2",
        "fwh3",
        "fwh4",
        "fwl1",
        "fwl2",
        "fwl3",
        "fwl4",
        "cdrh1",
        "cdrh2",
        "cdrh3",
        "cdrl1",
        "cdrl2",
        "cdrl3",
    ]
    _macro_regions = {
        "hframework": set(["fwh1", "fwh2", "fwh3", "fwh4"]),
        "hcdrs": set(["cdrh1", "cdrh2", "cdrh3"]),
        "lframework": set(["fwl1", "fwl2", "fwl3", "fwl4"]),
        "lcdrs": set(["cdrl1", "cdrl2", "cdrl3"]),
    }
    _macro_regions.update(
        {
            "framework": _macro_regions["hframework"] | _macro_regions["lframework"],
            "cdrs": _macro_regions["hcdrs"] | _macro_regions["lcdrs"],
            "vh": _macro_regions["hcdrs"] | _macro_regions["hframework"],
            "vl": _macro_regions["lcdrs"] | _macro_regions["lframework"],
        }
    )

    _macro_regions.update({"fv": _macro_regions["vh"] | _macro_regions["vl"]})

    _macro_positions = {}

    def __init__(self, numbering_scheme="chothia", definition="chothia", NOT=False):
        self.NOT = NOT
        self.set_regions()
        self.positions = {"H": set(), "L": set()}
        self.numbering_scheme = numbering_scheme
        self.definition = definition
        self.exclude = {"H": set(), "L": set()}

    def set_regions(self, regions=[]):
        """
        Set the regions to be used. Will clear anything added using add regions.
        """
        if self.NOT:
            self.regions = self._macro_regions["fv"]
        else:
            self.regions = set()
        self.add_regions(regions)

    def add_regions(self, regions):
        """
        Add regions to the selection.
        """
        for region in regions:
            region = region.lower()
            if region in self._defined_regions:
                if self.NOT:
                    self.regions = self.regions - set([region])
                else:
                    self.regions.add(region)
            elif region in self._macro_regions:
                if self.NOT:
                    self.regions = self.regions - self._macro_regions[region]
                else:
                    self.regions = self.regions | self._macro_regions[region]
            elif region in self._macro_positions:  # e.g. interface positions
                raise AssertionError
            else:
                raise AssertionError

    def add_positions(self, positions, chain):
        for position in positions:
            index, insertion = position
            self.positions[chain].add((index, insertion))

    def exclude_positions(self, positions, chain):
        for position in positions:
            index, insertion = position
            self.exclude[chain].add((index, insertion))

    def accept(self, position, chain):
        if position in self.exclude[chain]:
            return
        if (
            get_region(position, chain, self.numbering_scheme, self.definition)
            in self.regions
            or position in self.positions[chain]
        ):
            return 1


# def _generate_global_code():
#    from anarci.schemes import number_chothia_heavy, number_chothia_light,number_kabat_heavy,number_kabat_light,number_martin_heavy, number_martin_light,number_imgt
#    import textwrap
#    index_to_imgt_state = {}
#    state_vector   = [ ((i,"m"), i-1) for i in xrange( 1, 129 ) ]
#    dummy_sequence = "-"*128
#    scheme_to_number = { ("chothia", "H"): number_chothia_heavy,
#                         ("chothia", "L"): number_chothia_light,
#                         ("kabat", "H"):   number_kabat_heavy,
#                         ("kabat", "L"):   number_kabat_light,
#                         ("martin", "H"):  number_martin_heavy,
#                         ("martin", "L"):  number_martin_light,
#                         ("imgt", "H"):    number_imgt,
#                         ("imgt", "L"):    number_imgt
#                        }
#    for scheme, chain in scheme_to_number:
#        numbering = scheme_to_number[(scheme, chain)](state_vector, dummy_sequence)
#        index_to_imgt_state[(scheme, chain)] = dict( (numbering[0][i][0][0], i) for i in xrange( len( numbering[0] ) ) )
#    print "_index_to_imgt_state = ", "\n".join(textwrap.wrap( str( index_to_imgt_state ), 150 ))
#    print "_regions = {'imgt':{}}"
#    print "_regions['imgt']['L'] = _regions['imgt']['H'] = '11111111111111111111111111222222222222333333333333333334444444444555555555555555555555555555555555555555666666666666677777777777'"
#    for definition in [ "kabat","chothia","contact","north" ]:
#        print "_regions['%s'] = {}"%definition
#        Index = index_to_imgt_state[("chothia", "L")]
#        L1 = regions_in_chothia[ "L1" ][ definition ]
#        lrstr =  "1"*(Index[L1[0]])
#        lrstr += "2"*(Index[L1[1]+1]-Index[L1[0]])
#        L2 = regions_in_chothia[ "L2" ][ definition ]
#        lrstr +=  "3"*(Index[L2[0]]-Index[L1[1]]-1)
#        lrstr += "4"*(Index[L2[1]+1]-Index[L2[0]])
#        L3 = regions_in_chothia[ "L3" ][ definition ]
#        lrstr +=  "5"*(Index[L3[0]]-Index[L2[1]]-1)
#        lrstr += "6"*(Index[L3[1]+1]-Index[L3[0]])
#        lrstr += "7"*( 128 - Index[L3[1]]-1 )
#        print "_regions['%s']['L'] = "%definition, "'%s'"%lrstr
#        Index = index_to_imgt_state[("chothia", "H")]
#        H1 = regions_in_chothia[ "H1" ][ definition ]
#        hrstr =  "1"*(Index[H1[0]])
#        hrstr += "2"*(Index[H1[1]+1]-Index[H1[0]])
#        H2 = regions_in_chothia[ "H2" ][ definition ]
#        hrstr +=  "3"*(Index[H2[0]]-Index[H1[1]]-1)
#        hrstr += "4"*(Index[H2[1]+1]-Index[H2[0]])
#        H3 = regions_in_chothia[ "H3" ][ definition ]
#        hrstr +=  "5"*(Index[H3[0]]-Index[H2[1]]-1)
#        hrstr += "6"*(Index[H3[1]+1]-Index[H3[0]])
#        hrstr += "7"*( 128 - Index[H3[1]]-1 )
#        print "_regions['%s']['H'] = "%definition, "'%s'"%hrstr
