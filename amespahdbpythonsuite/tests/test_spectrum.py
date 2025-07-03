#!/usr/bin/env python3
"""
test_spectrum.py

Test the spectrum.py module.
"""

from os.path import exists

import astropy.units as u
import importlib_resources
import matplotlib.pyplot as plt
import numpy as np
import pytest
from astropy.io import ascii

from amespahdbpythonsuite import mcfitted, observation, spectrum
from amespahdbpythonsuite.amespahdb import AmesPAHdb


@pytest.fixture(scope="module")
def test_transitions():
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
    transitions = db.gettransitionsbyuid([18, 73, 726, 2054, 223])
    transitions.cascade(6 * 1.603e-12, multiprocessing=False)
    transitions.shift(-15.0)
    return transitions


@pytest.fixture(scope="module")
def test_spectrum(test_transitions):
    return test_transitions.convolve(fwhm=15.0)


@pytest.fixture(scope="module")
def test_observations():
    file = (
        importlib_resources.files("amespahdbpythonsuite") / "resources/galaxy_spec.ipac"
    )
    obs = observation.Observation(file)
    obs.abscissaunitsto("1/cm")
    return obs


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
        file = (
            importlib_resources.files("amespahdbpythonsuite")
            / "resources/galaxy_spec.ipac"
        )
        tbl = ascii.read(file)
        spectrum = test_transitions.convolve(
            grid=1e4 / tbl["wavelength"],
            fwhm=15.0,
            gaussian=True,
            multiprocessing=False,
        )
        fit = spectrum.fit(tbl["flux"], tbl["flux_uncertainty"])
        assert fit.getmethod() == "NNLC"

    def test_fit_without_errors(self, test_transitions):
        file = (
            importlib_resources.files("amespahdbpythonsuite")
            / "resources/sample_data_NGC7023.tbl"
        )
        tbl = ascii.read(file)
        spectrum = test_transitions.convolve(
            grid=tbl["WAVELENGTH"].to("1/cm", equivalencies=u.spectral()),
            fwhm=15.0,
            gaussian=True,
            multiprocessing=False,
        )
        fit = spectrum.fit(tbl["FLUX"])
        assert fit.getmethod() == "NNLS"

    def test_fit_with_obs_with_errors(self, test_observations, test_transitions):
        spectrum = test_transitions.convolve(
            grid=test_observations.spectrum.spectral_axis.value,
            fwhm=15.0,
            gaussian=True,
            multiprocessing=False,
        )
        fit = spectrum.fit(test_observations)
        assert fit.getmethod() == "NNLC"

    def test_fit_with_obs_without_errors(self, test_transitions):
        file = (
            importlib_resources.files("amespahdbpythonsuite")
            / "resources/sample_data_NGC7023.tbl"
        )
        obs = observation.Observation(file)
        obs.abscissaunitsto("1/cm")
        spectrum = test_transitions.convolve(
            grid=obs.spectrum.spectral_axis.value,
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

    def test_mcfit(self, test_observations, test_transitions):
        spectrum = test_transitions.convolve(
            grid=1e4 / test_observations.spectrum.spectral_axis.value,
            fwhm=15.0,
            gaussian=True,
            multiprocessing=False,
        )
        mcfit = spectrum.mcfit(test_observations, samples=10)
        assert isinstance(mcfit, mcfitted.MCFitted)
        assert len(mcfit.mcfits) == 10

    def test_mcfit_multiprocessing(self, test_observations, test_transitions):
        spectrum = test_transitions.convolve(
            grid=1e4 / test_observations.spectrum.spectral_axis.value,
            fwhm=15.0,
            gaussian=True,
            multiprocessing=False,
        )
        mcfit = spectrum.mcfit(test_observations, samples=10, multiprocessing=True)
        assert isinstance(mcfit, mcfitted.MCFitted)
        assert len(mcfit.mcfits) == 10
