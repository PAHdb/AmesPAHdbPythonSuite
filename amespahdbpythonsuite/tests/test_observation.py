#!/usr/bin/env python3
"""
test_obsservationervation.py

Test the observation.py module.
"""

import pytest
import numpy as np
import matplotlib.pyplot as plt

from pkg_resources import resource_filename


from amespahdbpythonsuite import observation


@pytest.fixture(scope="module")
def test_observation():
    file = resource_filename(
        "amespahdbpythonsuite", "resources/sample_data_NGC7023.tbl"
    )
    return observation.Observation(file)


class TestObservation:
    """
    Test Observation class.

    """

    def test_instance(self):
        assert isinstance(observation.Observation(), observation.Observation)

    def test_plot(self, monkeypatch, test_observation):
        monkeypatch.setattr(plt, "show", lambda: None)
        test_observation.plot(show=True)

    def test_read_fits(self):
        file = "resources/sample_data_NGC7023.fits"
        path = resource_filename("amespahdbpythonsuite", file)
        assert isinstance(observation.Observation(path), observation.Observation)

    def test_read_ipac(self):
        file = "resources/sample_data_NGC7023.tbl"
        path = resource_filename("amespahdbpythonsuite", file)
        assert isinstance(observation.Observation(path), observation.Observation)

    def test_file_not_exists(self):
        with pytest.raises(FileNotFoundError) as pytest_wrapped_e:
            observation.Observation("file_does_not_exist")
            assert pytest_wrapped_e.type == FileNotFoundError

    def test_file_malformed(self):
        file = "resources/sample_malformed.fits"
        path = resource_filename("amespahdbpythonsuite", file)
        with pytest.raises(OSError) as pytest_wrapped_e:
            observation.Observation(path)
            assert pytest_wrapped_e.type == OSError

    def test_rebin(self, test_observation):
        g1 = test_observation.getgrid()
        np.seterr(invalid="ignore")
        test_observation.rebin(0.25, uniform=True)
        np.seterr(invalid="warn")
        g2 = test_observation.getgrid()
        assert g2[0] == g1[0] and g2[-1] == g1[-1] and g2[1] - g2[0] == 0.25
        test_observation.rebin(200, resolution=True)
        g3 = test_observation.getgrid()
        assert (
            g3[0] == g2[0]
            and g3[-1] == g2[-1]
            and np.isclose(g3[0] / (g3[1] - g3[0]), 200)
        )

    def test_setgridrange(self, test_observation):
        test_observation.setgridrange(10.0, 12.0)
        g = test_observation.getgrid()
        assert g.min() >= 10.0 and g.max() <= 12.0

    def test_getset(self):
        file = "resources/sample_data_NGC7023.tbl"
        path = resource_filename("amespahdbpythonsuite", file)
        obs1 = observation.Observation(path)
        o1 = obs1.get()
        assert o1["type"] == "Observation"
        obs2 = observation.Observation()
        obs2.set(o1)
        o2 = obs2.get()
        assert o2["type"] == "Observation"

    def test_abscissaunitsto(self, test_observation):
        test_observation.abscissaunitsto("1/cm")
        assert test_observation.spectrum.spectral_axis.unit == "1/cm"
