#!/usr/bin/env python3
"""
test_spectrum.py

Test the spectrum.py module.
"""

import pytest
from os.path import exists
import numpy as np
import matplotlib.pyplot as plt
from pkg_resources import resource_filename

from astropy.io import ascii

from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite import spectrum, observation, mcfitted


@pytest.fixture(scope="module")
def test_transitions():
    xml = "resources/pahdb-theoretical_cutdown.xml"
    db = AmesPAHdb(
        filename=resource_filename("amespahdbpythonsuite", xml),
        check=False,
        cache=False,
        update=False,
    )
    transitions = db.gettransitionsbyuid([18, 73, 726, 2054, 223])
    transitions.cascade(6 * 1.603e-12, multiprocessing=False)
    transitions.shift(-15.0)
    return transitions


@pytest.fixture(scope="module")
def test_spectrum(test_transitions):
    return test_transitions.convolve(fwhm=15.0)


@pytest.fixture(scope="module")
def test_path(tmp_path_factory):
    d = tmp_path_factory.mktemp("test_spectrum")
    return f"{d}/result"


class TestSpectrum:
    """
    Test Spectrum class.

    """

    def test_instance(self):
        assert isinstance(spectrum.Spectrum(), spectrum.Spectrum)

    def test_plot(self, monkeypatch, test_spectrum):
        monkeypatch.setattr(plt, "show", lambda: None)
        test_spectrum.plot(show=True)

    def test_normalization(self, test_spectrum):
        test_spectrum.normalize()
        assert test_spectrum.data[73].max() == 1.0

    def test_fit_with_errors(self, test_transitions):
        file = resource_filename("amespahdbpythonsuite", "resources/galaxy_spec.ipac")
        f = ascii.read(file)
        spectrum = test_transitions.convolve(
            grid=[1e4 / x for x in f["wavelength"]],
            fwhm=15.0,
            gaussian=True,
            multiprocessing=False,
        )
        fit = spectrum.fit(f["flux"], f["flux_uncertainty"])
        assert fit.getmethod() == "NNLC"

    def test_fit_without_errors(self, test_transitions):
        file = resource_filename(
            "amespahdbpythonsuite", "resources/sample_data_NGC7023.tbl"
        )
        f = ascii.read(file)
        spectrum = test_transitions.convolve(
            grid=[1e4 / x for x in f["WAVELENGTH"]],
            fwhm=15.0,
            gaussian=True,
            multiprocessing=False,
        )
        fit = spectrum.fit(f["FLUX"])
        assert fit.getmethod() == "NNLS"

    def test_fit_with_obs_with_errors(self, test_transitions):
        file = resource_filename("amespahdbpythonsuite", "resources/galaxy_spec.ipac")
        obs = observation.Observation(file)
        spectrum = test_transitions.convolve(
            grid=1e4 / obs.spectrum.spectral_axis.value,
            fwhm=15.0,
            gaussian=True,
            multiprocessing=False,
        )
        fit = spectrum.fit(obs)
        assert fit.getmethod() == "NNLC"

    def test_fit_with_obs_without_errors(self, test_transitions):
        file = resource_filename(
            "amespahdbpythonsuite", "resources/sample_data_NGC7023.tbl"
        )
        obs = observation.Observation(file)
        spectrum = test_transitions.convolve(
            grid=1e4 / obs.spectrum.spectral_axis.value,
            fwhm=15.0,
            gaussian=True,
            multiprocessing=False,
        )
        fit = spectrum.fit(obs)
        assert fit.getmethod() == "NNLS"

    def test_resample(self, test_spectrum):
        g = np.arange(1000, 1600, 5)
        test_spectrum.resample(g)
        assert test_spectrum.grid.min() == g[0] and test_spectrum.grid.max() == g[-1]

    def test_getset(self, test_spectrum):
        s1 = test_spectrum.get()
        assert s1["type"] == "Spectrum"
        spec2 = spectrum.Spectrum()
        spec2.set(s1)
        s2 = spec2.get()
        assert s2["type"] == "Spectrum"

    def test_write_spectrum(self, test_spectrum, test_path):
        test_spectrum.write(f"{test_path}.tbl")
        assert exists(f"{test_path}.tbl")

    def test_mcfit(self, test_transitions):
        file = resource_filename("amespahdbpythonsuite", "resources/galaxy_spec.ipac")
        obs = observation.Observation(file)
        spectrum = test_transitions.convolve(
            grid=1e4 / obs.spectrum.spectral_axis.value,
            fwhm=15.0,
            gaussian=True,
            multiprocessing=False,
        )
        mcfit = spectrum.mcfit(obs, samples=10)
        assert isinstance(mcfit, mcfitted.MCfitted)
        assert len(mcfit.mcfits) == 10
