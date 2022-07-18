#!/usr/bin/env python3
"""
test_mcfitted.py

Test the mcfitted.py module.
"""

import pytest
from os.path import exists
import matplotlib.pyplot as plt

from pkg_resources import resource_filename

from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite import observation, mcfitted


@pytest.fixture(scope="module")
def test_mcfitted():
    xml = "resources/pahdb-theoretical_cutdown.xml"
    db = AmesPAHdb(
        filename=resource_filename("amespahdbpythonsuite", xml),
        check=False,
        cache=False,
        update=False,
    )
    uids = [18, 73, 726, 2054, 223]
    transitions = db.gettransitionsbyuid(uids)
    transitions.cascade(6 * 1.603e-12, multiprocessing=False)
    transitions.shift(-15.0)
    obs = observation.Observation(
        resource_filename("amespahdbpythonsuite", "resources/galaxy_spec.ipac")
    )
    spectrum = transitions.convolve(
        grid=1e4 / obs.getgrid(), fwhm=15.0, gaussian=True, multiprocessing=False
    )

    return spectrum.mcfit(obs, nsamples=10)


@pytest.fixture(scope="module")
def test_path(tmp_path_factory):
    d = tmp_path_factory.mktemp("test_mcfitted")
    print(d)
    return f"{d}/"


class TestMcfitted:
    """
    Test Spectrum class.

    """

    def test_instance(self):
        assert isinstance(mcfitted.MCfitted(), mcfitted.MCfitted)

    def test_plot(self, monkeypatch, test_mcfitted):
        monkeypatch.setattr(plt, "show", lambda: None)
        test_mcfitted.plot(show=True)

    def test_stats(self, test_mcfitted, test_path):
        test_mcfitted.getstats(save=True)
        assert exists('mcfitted_statistics.txt')

    def test_plot_charge(self, test_mcfitted, test_path):
        test_mcfitted.plot(
            wavelength=True,
            charge=True,
            save=test_path,
        )
        assert exists(f"{test_path}mc_charge_breakdown.pdf")

    def test_plot_size(self, test_mcfitted, test_path):
        test_mcfitted.plot(
            wavelength=True,
            size=True,
            save=test_path,
        )
        assert exists(f"{test_path}mc_size_breakdown.pdf")

    def test_plot_composition(self, test_mcfitted, test_path):
        test_mcfitted.plot(
            wavelength=True,
            composition=True,
            save=test_path,
        )
        assert exists(f"{test_path}mc_composition_breakdown.pdf")
