#!/usr/bin/env python3
"""
test_coadded.py

Test the coadded.py module.
"""

import pytest
import numpy as np
import matplotlib.pyplot as plt
from pkg_resources import resource_filename


from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite import coadded


@pytest.fixture(scope="module")
def test_coadded_result():
    file = "resources/coadded_test_data.npy"
    spec = np.load(resource_filename("amespahdbpythonsuite", file))
    return spec


@pytest.fixture(scope="module")
def test_coadded():
    xml = "resources/pahdb-theoretical_cutdown.xml"
    db = AmesPAHdb(
        filename=resource_filename("amespahdbpythonsuite", xml),
        check=False,
        cache=False,
        update=False,
    )
    trans = db.gettransitionsbyuid([18, 73])
    spec = trans.convolve(fwhm=15.0, gaussian=True)
    return spec.coadd(weights={18: 1.0, 73: 2.0}, average=True)


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
