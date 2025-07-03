#!/usr/bin/env python3
"""
test_coadded.py

Test the coadded.py module.
"""

from os.path import exists

import importlib_resources
import matplotlib.pyplot as plt
import numpy as np
import pytest

from amespahdbpythonsuite import coadded
from amespahdbpythonsuite.amespahdb import AmesPAHdb


@pytest.fixture(scope="module")
def test_coadded_result():
    spec = np.load(
        importlib_resources.files("amespahdbpythonsuite")
        / "resources/coadded_test_data.npy"
    )
    return spec


@pytest.fixture(scope="module")
def test_coadded():
    xml = (
        importlib_resources.files("amespahdbpythonsuite")
        / "resources/pahdb-theoretical_cutdown.xml"
    )
    db = AmesPAHdb(
        filename=xml,
        check=False,
        cache=False,
        update=False,
    )
    trans = db.gettransitionsbyuid([18, 73])
    spec = trans.convolve(fwhm=15.0, gaussian=True)
    return spec.coadd(weights={18: 1.0, 73: 2.0}, average=True)


@pytest.fixture(scope="module")
def test_path(tmp_path_factory):
    d = tmp_path_factory.mktemp("test_coadded")
    return f"{d}/result"


class TestCoadded:
    """
    Test Coadded class.

    """

    def test_instance(self, test_coadded, test_coadded_result):
        np.testing.assert_allclose(test_coadded_result, test_coadded.data[0])

    def test_plot(self, monkeypatch, test_coadded):
        monkeypatch.setattr(plt, "show", lambda: None)
        test_coadded.plot(show=True)

    def test_getset(self, test_coadded):
        c1 = test_coadded.get()
        assert c1["type"] == "Coadded"
        coadded2 = coadded.Coadded()
        coadded2.set(c1)
        c2 = coadded2.get()
        assert c2["type"] == "Coadded"

    def test_write_coadded(self, test_coadded, test_path):
        test_coadded.write(f"{test_path}.tbl")
        assert exists(f"{test_path}.tbl")
