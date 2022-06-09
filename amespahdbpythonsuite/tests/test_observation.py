#!/usr/bin/env python3
"""
test_observation.py

Test the observation.py module.
"""

import pytest

from pkg_resources import resource_filename


from amespahdbpythonsuite import observation


class TestObservation():
    """
    Test Observation class.

    """

    def test_read_fits(self):
        file = 'resources/sample_data_NGC7023.fits'
        path = resource_filename('amespahdbpythonsuite', file)
        assert isinstance(observation.Observation(path),
                          observation.Observation)

    def test_read_ipac(self):
        file = 'resources/sample_data_NGC7023.tbl'
        path = resource_filename('amespahdbpythonsuite', file)
        assert isinstance(observation.Observation(path),
                          observation.Observation)

    def test_file_not_exists(self):
        with pytest.raises(FileNotFoundError) as pytest_wrapped_e:
            observation.Observation('file_does_not_exist')
            assert pytest_wrapped_e.type == FileNotFoundError

    def test_file_malformed(self):
        file = 'resources/sample_malformed.fits'
        path = resource_filename('amespahdbpythonsuite', file)
        with pytest.raises(OSError) as pytest_wrapped_e:
            observation.Observation(path)
            assert pytest_wrapped_e.type == OSError

    def test_getset(self):
        file = 'resources/sample_data_NGC7023.tbl'
        path = resource_filename('amespahdbpythonsuite', file)
        obs1 = observation.Observation(path)
        o1 = obs1.get()
        assert(o1['type'] == 'Observation')
        obs2 = observation.Observation()
        obs2.set(o1)
        o2 = obs2.get()
        assert(o2['type'] == 'Observation')
