#!/usr/bin/env python3
"""
test_spectrum.py

Test the spectrum.py module.
"""

import pytest
import numpy as np
import matplotlib.pyplot as plt
from pkg_resources import resource_filename

from astropy.io import ascii

from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite import spectrum


@pytest.fixture(scope="module")
def pahdb_theoretical():
    xml = 'resources/pahdb-theoretical_cutdown.xml'
    pahdb = AmesPAHdb(filename=resource_filename('amespahdbpythonsuite', xml),
                      check=False, cache=False)
    return pahdb


@pytest.fixture(scope="module")
def test_spec():

    file = 'resources/coadded_test_data.npy'
    spec = np.load(resource_filename('amespahdbpythonsuite', file))

    return spec


class TestSpectrum():
    """
    Test Spectrum class.

    """

    def test_fit(self, pahdb_theoretical, monkeypatch, test_spec):
        # Read input spectrum.
        spec = resource_filename(
            'amespahdbpythonsuite', 'resources/galaxy_spec.ipac')
        f = ascii.read(spec)
        wave = f['wavelength']
        flux = f['flux']
        sigma = f['sigma']
        waven = [1e4 / x for x in wave]
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18, 73, 726, 2054, 223]
        # Retrieve the transitions from the database for the subset of PAHs.
        transitions = pahdb.gettransitionsbyuid(uids)
        # Set emission model.
        transitions.cascade(6 * 1.603e-12, multiprocessing=False)
        # Shift data 15 wavenumber to the red.
        transitions.shift(-15.0)
        # convolve the transitions into a spectrum.
        spectrum = transitions.convolve(
            grid=waven, fwhm=15.0, gaussian=True, multiprocessing=False)
        # fit the spectrum.
        fit = spectrum.fit(flux, sigma)
        # Assert results.
        assert spectrum.uids == uids
        assert fit.uids == [73, 2054, 223]
        np.testing.assert_array_almost_equal(fit.grid, np.array(waven))
        # Check plotting function.
        monkeypatch.setattr(plt, 'show', lambda: None)
        spectrum.plot()
        coadded = fit.coadd(weights=fit.weights, average=False)
        np.testing.assert_equal(test_spec, coadded['data'])

    def test_getset(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18, 73]
        trans = pahdb.gettransitionsbyuid(uids)
        spec1 = trans.convolve(fwhm=15.0)
        d1 = spec1.get()
        assert(d1['type'] == 'Spectrum')
        spec2 = spectrum.Spectrum()
        spec2.set(d1)
        d2 = spec2.get()
        assert(d2['type'] == 'Spectrum')
