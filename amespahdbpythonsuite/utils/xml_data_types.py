#!/usr/bin/env python3
"""
xml_data_types.py

Desired data types when extracting data from a PAHdb XML file.
"""

metadata_keys = {
    'references': str,
    'comments': str,
    'formula': str,
    'charge': int,
    'symmetry': str,
    'weight': float,
    'total_e': float,
    'vib_e': float,
    'method': str,
    'n_solo': int,
    'n_duo': int,
    'n_trio': int,
    'n_quartet': int,
    'n_quintet': int,
    'n_ch2': int,
    'n_chx': int,
}

geometry_keys = {
    'position': int,
    # 'uid': int,
    # 'type': int,
    'x': float,
    'y': float,
    'z': float,
}

transition_keys = {
    'frequency': float,
    # 'uid': int,
    'scale': float,
    'intensity': float,
    'symmetry': str,
}

laboratory_keys = {
    'frequency': float,
    'intensity': float,
}
