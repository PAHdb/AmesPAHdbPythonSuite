#!/usr/bin/env python3
"""
test_geometry.py

Test the geometry.py module.
"""

import pytest
from pkg_resources import resource_filename

import numpy as np

from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite import geometry


@pytest.fixture(scope="module")
def geometry_test():
    xml = "resources/pahdb-theoretical_cutdown.xml"
    db = AmesPAHdb(
        filename=resource_filename("amespahdbpythonsuite", xml),
        check=False,
        cache=False,
        update=False,
    )
    g = db.getgeometrybyuid([18, 73, 726, 2054, 223])
    return g


@pytest.fixture(scope="module")
def test_masses():
    file = "resources/test_geometry_masses.npy"
    return np.load(resource_filename("amespahdbpythonsuite", file), allow_pickle=True)


@pytest.fixture(scope="module")
def test_diagonalized():
    file = "resources/test_geometry_diagonalized.npy"
    return np.load(resource_filename("amespahdbpythonsuite", file))


@pytest.fixture(scope="module")
def test_tensor():
    file = "resources/test_geometry_tensor.npy"

    return np.load(resource_filename("amespahdbpythonsuite", file))


@pytest.fixture(scope="module")
def test_nrings():
    file = "resources/test_geometry_rings.npy"
    return np.load(resource_filename("amespahdbpythonsuite", file), allow_pickle=True)


@pytest.fixture(scope="module")
def test_areas():
    file = "resources/test_geometry_areas.npy"
    return np.load(resource_filename("amespahdbpythonsuite", file), allow_pickle=True)


class TestGeometry:
    """
    Test Geometry class.

    """

    def test_instance(self):
        assert isinstance(geometry.Geometry(), geometry.Geometry)

    def test_mass(self, geometry_test, test_masses):
        assert geometry_test.mass() == test_masses

    def test_inertia(self, geometry_test, test_tensor):
        inertia = geometry_test.inertia()[18]
        np.testing.assert_allclose(inertia, test_tensor)

    def test_diagonalize(self, geometry_test, test_diagonalized):
        geometry_test.diagonalize()
        g = geometry_test.get()
        x = [d["x"] for d in g["data"][18]]
        np.testing.assert_allclose(x, test_diagonalized)

    def test_rings(self, geometry_test, test_nrings):
        rings = geometry_test.rings()
        assert rings[726] == test_nrings

    def test_area(self, geometry_test, test_areas):
        areas = geometry_test.area()
        assert areas == test_areas

    def test_getset(self, geometry_test):
        g1 = geometry_test.get()
        assert g1["type"] == "Geometry"
        geometry2 = geometry.Geometry()
        geometry2.set(g1)
        g2 = geometry2.get()
        assert g2["type"] == "Geometry"
