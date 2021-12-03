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
from amespahdbpythonsuite import transitions


@pytest.fixture(scope="module")
def pahdb_theoretical():
    xml = 'resources/pahdb-theoretical_cutdown.xml'
    pahdb = AmesPAHdb(filename=resource_filename('amespahdbpythonsuite', xml),
                      check=False, cache=False)
    return pahdb


@pytest.fixture(scope="module")
def test_spec():
    tspec1 = np.array([0.0, 0.0,
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

    tspec2 = np.array([1.67804028e-16, 2.81683088e-16,
                       6.20993619e-16, 2.66015181e-15,
                       2.02113387e-15, 1.55593909e-15,
                       3.24139847e-15, 7.93277653e-15,
                       2.92710663e-15, 2.39183719e-15,
                       3.26750798e-15, 1.80596235e-15,
                       1.67358352e-15, 2.21047413e-15,
                       3.60529180e-15, 7.24592115e-15,
                       1.86989149e-14, 3.80075950e-14,
                       1.91392555e-14, 8.74063194e-15,
                       6.35264569e-15, 4.70897652e-15,
                       2.68812711e-15, 1.76447652e-15,
                       1.31117613e-15, 1.06195248e-15,
                       9.23434317e-16, 8.60352326e-16,
                       8.63953102e-16, 9.46000302e-16,
                       1.14762718e-15, 1.57221725e-15,
                       2.47743532e-15, 4.43921295e-15,
                       7.41744509e-15, 6.92983696e-15,
                       4.01907480e-15, 2.30061372e-15])

    tspec3 = np.array([2.07601509e-17, 3.69268953e-17,
                       9.08725641e-17, 1.08806954e-15,
                       4.68912433e-16, 3.13743924e-16,
                       1.12159613e-15, 2.24033681e-14,
                       7.58990395e-16, 6.80496592e-16,
                       6.04732430e-15, 4.53334019e-16,
                       3.89219648e-16, 5.35198632e-16,
                       9.37489208e-16, 2.16206714e-15,
                       8.76050847e-15, 1.41669509e-13,
                       9.58459994e-15, 3.19883865e-15,
                       4.19385441e-15, 3.36560048e-15,
                       9.98928375e-16, 6.08179823e-16,
                       4.53272698e-16, 3.74592550e-16,
                       3.34212960e-16, 3.20319997e-16,
                       3.31712348e-16, 3.76324233e-16,
                       4.78183370e-16, 7.04403991e-16,
                       1.27794166e-15, 3.23808693e-15,
                       1.32345274e-14, 1.02742723e-14,
                       2.91666426e-15, 1.31679190e-15])

    return tspec1, tspec2, tspec3


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

    def test_set(self, capsys, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        # Test transitions set and data class set methods.
        transitions.Transitions(type='invalid',
                                version=pahdb._AmesPAHdb__data['version'],
                                pahdb=pahdb._AmesPAHdb__data)
        captured = capsys.readouterr()
        assert 'DATABASE MISMATCH' in captured.out

        transitions.Transitions(type='theoretical',
                                version='9999',
                                pahdb=pahdb._AmesPAHdb__data)
        captured = capsys.readouterr()
        assert 'VERSION MISMATCH' in captured.out

    def test_get(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18, 73]
        # Test transitions set and data class set methods.
        trans = transitions.Transitions(type='theoretical',
                                        version=pahdb._AmesPAHdb__data['version'],
                                        uids=uids,
                                        data=pahdb._AmesPAHdb__getkeybyuids('transitions', uids),
                                        pahdb=pahdb._AmesPAHdb__data,
                                        model={'type': 'zerokelvin_m',
                                               'temperature': 0.0,
                                               'description': ''},
                                        units={'abscissa': {'unit': 1,
                                                            'str': 'frequency [wavenumber]'},
                                               'ordinate': {'unit': 2,
                                                            'str': 'integrated cross-section' + '[km/mol]'}})

        test_dict = trans.get()
        key_list = ['type', 'database', 'version', 'data', 'uids', 'model', 'units', 'shift']
        assert list(test_dict.keys()) == key_list
        assert test_dict['type'] == 'Transitions'
        assert test_dict['database'] == 'theoretical'
        assert test_dict['version'] == pahdb._AmesPAHdb__data['version']

    def test_intersect(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18, 73, 726, 2054, 223]
        sub_uids = [18, 223]
        trans = pahdb.gettransitionsbyuid(uids)
        trans.intersect(sub_uids)
        assert list(trans.uids) == sub_uids
        assert list(trans.data.keys()) == sub_uids

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
        # Apply cascade emission model at 6 eV.
        trans.cascade(6 * 1.603e-12, multiprocessing=False)
        # Assert attained Tmax.
        dtest_a = next((sub for sub in trans.model['temperatures'] if sub['uid'] == uids[0]), None)
        # Retrieve dictionary at 3068.821 frequency.
        dtest_b = [x for x in trans.data[18] if x['frequency'] == 3068.821][0]
        # Assert attained intensity and temperature.
        assert dtest_a['temperature'] == 1279.7835033561428
        assert dtest_b['intensity'] == 1.6710637100014386e-12

    def test_partial_cascade(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18]
        trans_multi = pahdb.gettransitionsbyuid(uids)
        # Assert Tmax and intensity at 3068.821 frequency.
        intf, tmax = trans_multi._cascade_em_model(6 * 1.603e-12, uids[0])
        test_i = [x for x in intf[uids[0]] if x['frequency'] == 3068.821][0]
        assert tmax['temperature'] == 1279.7835033561428
        assert test_i['intensity'] == 1.6710637100014386e-12

    def test_convolve(self, pahdb_theoretical, test_spec):
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18]
        trans = pahdb.gettransitionsbyuid(uids)
        # Create wavenumber grid.
        waven = [1e4 / x for x in np.arange(5, 20, 0.4)]
        # Apply cascade emission model at 6 eV.
        trans.cascade(6 * 1.603e-12, multiprocessing=False)
        # Apply shift.
        trans.shift(-15.0)
        # Get spectum object
        spec1 = trans.convolve(grid=waven, fwhm=15.0, gaussian=True, multiprocessing=False)
        spec2 = trans.convolve(grid=waven, fwhm=15.0, drude=True, multiprocessing=False)
        spec3 = trans.convolve(grid=waven, fwhm=15.0, multiprocessing=False)
        # Expected gaussian convolved transitions.
        assert spec1.model['temperatures'][0]['temperature'] == 1279.7835033561428
        np.testing.assert_array_almost_equal(test_spec[0], spec1.data[18])
        np.testing.assert_array_almost_equal(test_spec[1], spec2.data[18])
        np.testing.assert_array_almost_equal(test_spec[2], spec3.data[18])

    def test_partial_convolve(self, pahdb_theoretical, test_spec):
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18]
        trans_multi = pahdb.gettransitionsbyuid(uids)
        # Create wavenumber grid.
        waven = [1e4 / x for x in np.arange(5, 20, 0.4)]
        # Apply cascade emission model at 6 eV.
        trans_multi.cascade(6 * 1.603e-12, multiprocessing=True)
        # Apply shift.
        trans_multi.shift(-15.0)
        # Get spectum object
        spec = trans_multi._get_intensities(npoints=len(waven),
                                            xmin=np.min(waven),
                                            xmax=np.max(waven),
                                            clip=3,
                                            width=0.5 * 15.0 / np.sqrt(2.0 * np.log(2.0)),
                                            x=np.asarray(waven),
                                            gaussian=True,
                                            drude=False,
                                            uid=uids[0])
        np.testing.assert_array_almost_equal(test_spec[0], spec[18])

    def test_plot_transitions(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18]
        trans = pahdb.gettransitionsbyuid(uids)
        # Call plot method.
        trans.plot()
