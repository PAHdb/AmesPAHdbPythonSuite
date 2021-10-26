#!/usr/bin/env python3
"""
test_amespahdb.py

Test the Transitions.py module.
"""

import amespahdbpythonsuite

from amespahdbpythonsuite.geometry import Geometry


class TestTransitions():
    """
    Test AmesPAHdb class.

    """
    def test_instance(self):
        geo = Geometry()
        assert isinstance(geo, amespahdbpythonsuite.geometry.Geometry)
