#!/usr/bin/env python3
"""
test_geometry.py

Test the geometry.py module.
"""

import pytest
from pkg_resources import resource_filename
from os.path import exists
import matplotlib.pyplot as plt
import numpy as np

from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite import geometry


@pytest.fixture(scope="module")
def test_geometry():
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
def test_path(tmp_path_factory):
    d = tmp_path_factory.mktemp("test_geometry")
    return f"{d}/result"


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

    def test_plot(self, monkeypatch, test_geometry):
        monkeypatch.setattr(plt, "show", lambda: None)
        test_geometry.plot(18, show=True)

    def test_structure(self, test_geometry):
        test_geometry.structure(18, show=False)

    def test_structure_save(self, test_path, test_geometry):
        test_geometry.structure(18, transparent=True, save=test_path)
        assert exists(f"{test_path}.png")

    def test_mass(self, test_geometry, test_masses):
        assert test_geometry.mass() == test_masses

    def test_inertia(self, test_geometry, test_tensor):
        inertia = test_geometry.inertia()[18]
        np.testing.assert_allclose(inertia, test_tensor)

    def test_diagonalize(self, test_geometry, test_diagonalized):
        test_geometry.diagonalize()
        g = test_geometry.get()
        x = [d["x"] for d in g["data"][18]]
        np.testing.assert_allclose(x, test_diagonalized)

    def test_rings(self, test_geometry, test_nrings):
        rings = test_geometry.rings()
        assert rings[726] == test_nrings

    def test_area(self, test_geometry, test_areas):
        areas = test_geometry.area()
        assert areas == test_areas

    def test_bec(self, test_geometry):
        assert test_geometry.bec()[18] == "333333"

    def test_getset(self, test_geometry):
        g1 = test_geometry.get()
        assert g1["type"] == "Geometry"
        geometry2 = geometry.Geometry()
        geometry2.set(g1)
        g2 = geometry2.get()
        assert g2["type"] == "Geometry"
