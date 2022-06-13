#!/usr/bin/env python3
"""
test_coadded.py

Test the coadded.py module.
"""

import pytest
import numpy as np
from pkg_resources import resource_filename


from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite import coadded


@pytest.fixture(scope="module")
def pahdb_theoretical():
    xml = "resources/pahdb-theoretical_cutdown.xml"
    pahdb = AmesPAHdb(
        filename=resource_filename("amespahdbpythonsuite", xml),
        check=False,
        cache=False,
        update=False,
    )
    return pahdb


@pytest.fixture(scope="module")
def test_coadded():

    file = "resources/coadded_test_data.npy"
    spec = np.load(resource_filename("amespahdbpythonsuite", file))

    return spec


class TestCoadded:
    """
    Test Coadded class.

    """

    def test_initialization(self, pahdb_theoretical, test_coadded):
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18, 73]
        # Retrieve the transitions from the database for the subset of PAHs.
        transitions = pahdb.gettransitionsbyuid(uids)
        # convolve the transitions into a spectrum.
        spectrum = transitions.convolve(fwhm=15.0, gaussian=True, multiprocessing=False)
        coadded = spectrum.coadd(weights={18: 1.0, 73: 2.0}, average=True)
        np.testing.assert_allclose(test_coadded, coadded.data[0])

    def test_getset(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18, 73]
        trans = pahdb.gettransitionsbyuid(uids)
        spec = trans.convolve(fwhm=15.0)
        coadded1 = spec.coadd(weights={18: 1.0, 73: 2.0})
        c1 = coadded1.get()
        assert c1["type"] == "Coadded"
        coadded2 = coadded.Coadded()
        coadded2.set(c1)
        c2 = coadded2.get()
        assert c2["type"] == "Coadded"
