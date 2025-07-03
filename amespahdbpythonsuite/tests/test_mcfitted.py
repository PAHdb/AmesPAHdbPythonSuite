#!/usr/bin/env python3
"""
test_mcfitted.py

Test the mcfitted.py module.
"""

from os.path import exists

import importlib_resources
import matplotlib.pyplot as plt
import pytest

from amespahdbpythonsuite import mcfitted, observation
from amespahdbpythonsuite.amespahdb import AmesPAHdb


@pytest.fixture(scope="module")
def test_mcfitted():
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
    uids = [18, 73, 726, 2054, 223]
    transitions = db.gettransitionsbyuid(uids)
    obs = observation.Observation(
        importlib_resources.files("amespahdbpythonsuite") / "resources/galaxy_spec.ipac"
    )
    obs.abscissaunitsto("1/cm")
    spectrum = transitions.convolve(
        grid=obs.getgrid(), fwhm=15.0, gaussian=True, multiprocessing=False
    )
    return spectrum.mcfit(obs, samples=10)


@pytest.fixture(scope="module")
def test_path(tmp_path_factory):
    d = tmp_path_factory.mktemp("test_mcfitted")
    return f"{d}/"


class TestMCFitted:
    """
    Test Spectrum class.

    """

    def test_instance(self):
        assert isinstance(mcfitted.MCFitted(), mcfitted.MCFitted)

    def test_getfit(self, test_mcfitted):
        fit = test_mcfitted.getfit()
        assert len(fit.keys()) == 4

    def test_getclasses(self, test_mcfitted):
        classes = test_mcfitted.getclasses()
        assert len(classes.keys()) == 8

    def test_getbreakdown(self, test_mcfitted):
        breakdown = test_mcfitted.getbreakdown()
        assert len(breakdown.keys()) == 14

    def test_plot(self, monkeypatch, test_mcfitted):
        monkeypatch.setattr(plt, "show", lambda: None)
        test_mcfitted.plot(show=True)

    def test_plot_charge(self, test_mcfitted, test_path):
        test_mcfitted.plot(
            wavelength=True, charge=True, save=True, output=test_path, ftype="pdf"
        )
        assert exists(f"{test_path}mc_charge_breakdown.pdf")

    def test_plot_size(self, test_mcfitted, test_path):
        test_mcfitted.plot(
            wavelength=True, size=True, save=True, output=test_path, ftype="pdf"
        )
        assert exists(f"{test_path}mc_size_breakdown.pdf")

    def test_plot_composition(self, test_mcfitted, test_path):
        test_mcfitted.plot(
            wavelength=True, composition=True, save=True, output=test_path, ftype="pdf"
        )
        assert exists(f"{test_path}mc_composition_breakdown.pdf")

    def test_write(self, test_path, test_mcfitted):
        test_mcfitted.write(filename=f"{test_path}mc_breakdown.tbl")
        assert exists(f"{test_path}mc_breakdown.tbl")
