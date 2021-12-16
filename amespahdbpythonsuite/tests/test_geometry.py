#!/usr/bin/env python3
"""
test_geometry.py

Test the geometry.py module.
"""

import amespahdbpythonsuite

from amespahdbpythonsuite.geometry import Geometry


class TestGeometry():
    """
    Test Geometry class.

    """
    def test_instance(self):
        geo = Geometry()
        assert isinstance(geo, amespahdbpythonsuite.geometry.Geometry)
