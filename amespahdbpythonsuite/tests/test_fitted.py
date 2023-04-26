#!/usr/bin/env python3
"""
test_fitted.py

Test the fitted.py module.
"""

import pytest
from os.path import exists
from pkg_resources import resource_filename

import matplotlib.pyplot as plt

from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite import observation, fitted


@pytest.fixture(scope="module")
def test_fitted():
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
    obs.abscissaunitsto("1/cm")
    spectrum = transitions.convolve(
        grid=obs.getgrid(), fwhm=15.0, gaussian=True, multiprocessing=False
    )

    return spectrum.fit(obs)


@pytest.fixture(scope="module")
def test_path(tmp_path_factory):
    d = tmp_path_factory.mktemp("test_fitted")
    return f"{d}/result"


class TestFitted:
    """
    Test Fitted class.

    """

    def test_instance(self):
        assert isinstance(fitted.Fitted(), fitted.Fitted)

    def test_method(self, test_fitted):
        assert test_fitted.getmethod() == "NNLC"

    def test_write_fitted(self, test_fitted, test_path):
        test_fitted.write(f"{test_path}.tbl")
        assert exists(f"{test_path}.tbl")

    def test_breakdown(self, test_fitted):
        assert list(test_fitted.getbreakdown().keys()) == [
            "solo",
            "duo",
            "trio",
            "quartet",
            "quintet",
            "anion",
            "neutral",
            "cation",
            "small",
            "large",
            "nitrogen",
            "pure",
            "nc",
        ]

    def test_error(self, test_fitted):
        assert list(test_fitted.geterror().keys()) == [
            "err",
            "e127",
            "e112",
            "e77",
            "e62",
            "e33",
        ]

    def test_plot(self, test_fitted, test_path):
        test_fitted.plot(
            wavelength=True,
            sigma=test_fitted.observation.uncertainty.array,
            save=True,
            output=test_path,
            ptype="UIDs",
            ftype="pdf",
            units=[
                test_fitted.observation.spectral_axis.unit.to_string(),
                test_fitted.observation.flux.unit.to_string(),
            ],
        )
        assert exists(f"{test_path}_UIDs.pdf")

    def test_plot_residual(self, test_fitted, test_path):
        test_fitted.plot(
            wavelength=True,
            residual=True,
            sigma=test_fitted.observation.uncertainty.array,
            save=True,
            output=test_path,
            ptype="residual",
            ftype="pdf",
            units=[
                test_fitted.observation.spectral_axis.unit.to_string(),
                test_fitted.observation.flux.unit.to_string(),
            ],
        )
        assert exists(f"{test_path}_residual.pdf")

    def test_plot_size(self, test_fitted, test_path):
        test_fitted.plot(
            wavelength=True,
            size=True,
            sigma=test_fitted.observation.uncertainty.array,
            save=True,
            output=test_path,
            ptype="size",
            ftype="pdf",
            units=[
                test_fitted.observation.spectral_axis.unit.to_string(),
                test_fitted.observation.flux.unit.to_string(),
            ],
        )
        assert exists(f"{test_path}_size.pdf")

    def test_plot_charge(self, test_fitted, test_path):
        test_fitted.plot(
            wavelength=True,
            charge=True,
            sigma=test_fitted.observation.uncertainty.array,
            save=True,
            output=test_path,
            ptype="charge",
            ftype="pdf",
            units=[
                test_fitted.observation.spectral_axis.unit.to_string(),
                test_fitted.observation.flux.unit.to_string(),
            ],
        )
        assert exists(f"{test_path}_charge.pdf")

    def test_plot_composition(self, test_fitted, test_path):
        test_fitted.plot(
            wavelength=True,
            composition=True,
            sigma=test_fitted.observation.uncertainty.array,
            save=True,
            output=test_path,
            ptype="composition",
            ftype="pdf",
            units=[
                test_fitted.observation.spectral_axis.unit.to_string(),
                test_fitted.observation.flux.unit.to_string(),
            ],
        )
        assert exists(f"{test_path}_composition.pdf")

    def test_plot_sizedistribution(self, test_fitted, monkeypatch):
        monkeypatch.setattr(plt, "show", lambda: None)
        test_fitted.plot(sizedistribution=True, show=True)

    def test_sizedistribution(self, test_fitted):
        h, edges = test_fitted.getsizedistribution()
        assert len(h) == 3 and len(edges) == 4

    def test_sort(self, test_fitted):
        assert test_fitted.sort() != test_fitted.sort(flux=True)

    def test_getset(self, test_fitted):
        f1 = test_fitted.get()
        assert f1["type"] == "Fitted"
        fit = fitted.Fitted()
        fit.set(f1)
        f2 = fit.get()
        assert f2["type"] == "Fitted"
