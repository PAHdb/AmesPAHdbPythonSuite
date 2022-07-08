#!/usr/bin/env python3
"""
test_mcfitted.py

Test the mcfitted.py module.
"""

import pytest
from os.path import exists
import numpy as np
import matplotlib.pyplot as plt

from pkg_resources import resource_filename

from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite import mcfitted, observation


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


class TestMcfitted:
    """
    Test Spectrum class.

    """

    def test_instance(self):
        assert isinstance(mcfitted.Mcfitted(), mcfitted.Mcfitted)
