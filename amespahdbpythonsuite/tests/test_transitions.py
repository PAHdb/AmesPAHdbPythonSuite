#!/usr/bin/env python3
"""
test_amespahdb.py

Test the Transitions.py module.
"""

import numpy as np
import pytest
from pkg_resources import resource_filename
import amespahdbpythonsuite

from amespahdbpythonsuite.amespahdb import AmesPAHdb


@pytest.fixture(scope="module")
def pahdb_theoretical():
    xml = 'resources/pahdb-theoretical_cutdown.xml'
    pahdb = AmesPAHdb(filename=resource_filename('amespahdbpythonsuite', xml),
                      check=False, cache=False)
    return pahdb


class TestTransitions():
    """
    Test AmesPAHdb class.

    """
    def test_instance(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18, 73, 726, 2054, 223]
        trans = pahdb.gettransitionsbyuid(uids)
        assert isinstance(trans, amespahdbpythonsuite.transitions.Transitions)

    def test_shift(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18]
        trans = pahdb.gettransitionsbyuid(uids)
        # Retrieve dictionary at 3068.821 frequency.
        dtest = [x for x in trans.data[18] if x['frequency'] == 3068.821][0]
        # Apply shift.
        trans.shift(-15.0)
        # Assert attained intensity.
        assert dtest['frequency'] == 3053.821

    def test_fixedtemperature(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18]
        trans = pahdb.gettransitionsbyuid(uids)
        # Apply fixed temperature emission model at 600 K.
        trans.fixedtemperature(600)
        # Retrieve dictionary at 3068.821 frequency.
        dtest = [x for x in trans.data[18] if x['frequency'] == 3068.821][0]
        # Assert attained intensity.
        assert dtest['intensity'] == 6.420001406551514e-14

    def test_calculatedtemperature(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18]
        trans = pahdb.gettransitionsbyuid(uids)
        # Apply calculatedtemperature emission model at 6 eV.
        trans.calculatedtemperature(6 * 1.603e-12)
        # Assert attained Tmax.
        dtest = next((sub for sub in trans.model['temperatures'] if sub['uid'] == uids[0]), None)
        assert dtest['temperature'] == 1279.7835033561428

    def test_cascade(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18]
        trans = pahdb.gettransitionsbyuid(uids)
        trans_multi = pahdb.gettransitionsbyuid(uids)
        # Apply cascade emission model at 6 eV.
        trans.cascade(6 * 1.603e-12, multiprocessing=False)
        trans_multi.cascade(6 * 1.603e-12, multiprocessing=True)
        # Assert attained Tmax.
        dtest_a = next((sub for sub in trans.model['temperatures'] if sub['uid'] == uids[0]), None)
        dtest_am = next((sub for sub in trans_multi.model['temperatures'] if sub['uid'] == uids[0]), None)
        # Retrieve dictionary at 3068.821 frequency.
        dtest_b = [x for x in trans.data[18] if x['frequency'] == 3068.821][0]
        dtest_bm = [x for x in trans_multi.data[18] if x['frequency'] == 3068.821][0]
        # Assert attained intensity and temperature.
        assert dtest_a['temperature'] == 1279.7835033561428
        assert dtest_a['temperature'] == dtest_am['temperature']
        assert dtest_b['intensity'] == 1.6710637100014386e-12
        assert dtest_b['intensity'] == dtest_bm['intensity']

    def test_convolve(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18]
        trans = pahdb.gettransitionsbyuid(uids)
        trans_multi = pahdb.gettransitionsbyuid(uids)
        # Create wavenumber grid.
        waven = [1e4 / x for x in np.arange(5, 20, 0.4)]
        # Apply cascade emission model at 6 eV.
        trans.cascade(6 * 1.603e-12, multiprocessing=False)
        trans_multi.cascade(6 * 1.603e-12, multiprocessing=True)
        # Apply shift.
        trans.shift(-15.0)
        trans_multi.shift(-15.0)
        # Get spectum object
        spec = trans.convolve(grid=waven, fwhm=15.0, gaussian=True, multiprocessing=False)
        spec_multi = trans_multi.convolve(grid=waven, fwhm=15.0, gaussian=True, multiprocessing=True)
        # Expected gaussian convolved transitions.
        tspec = np.array([0.0, 0.0,
                          3.4411188697810187e-128, 1.2190842078939826e-20,
                          8.919790520493667e-26, 1.2275729151921775e-24,
                          3.4503464412612805e-16, 3.3255951213183884e-14,
                          2.2845397886608066e-19, 7.376715353567142e-20,
                          8.638506960894798e-15, 1.690177885781663e-30,
                          5.302322889355252e-66, 1.828462880535185e-107,
                          2.7129757078463156e-62, 1.6960157922699966e-33,
                          5.05042952727869e-18, 2.090532525483598e-13,
                          1.2056029704312457e-17, 1.5888762280459496e-16,
                          2.964176543593395e-15, 2.3015310890351513e-15,
                          1.2528641227452116e-20, 3.44802322919447e-30,
                          2.5192831665285025e-43, 1.9694271817882437e-59,
                          5.331856658440902e-78, 9.565651820272891e-68,
                          2.452346086154861e-52, 4.2447436371370704e-40,
                          1.2209172829656069e-30, 1.2771573679264533e-23,
                          9.617501977937132e-19, 9.467817572725245e-16,
                          2.0545192799741963e-14, 1.5543597485079575e-14,
                          6.135216019051803e-16, 1.8019125715086965e-18])
        assert spec.model['temperatures'][0]['temperature'] == 1279.7835033561428
        assert spec.model['temperatures'][0]['temperature'] == spec_multi.model['temperatures'][0]['temperature']
        np.testing.assert_array_almost_equal(tspec, spec.data[18])
        np.testing.assert_array_almost_equal(tspec, spec_multi.data[18])
