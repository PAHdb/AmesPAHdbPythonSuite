#!/usr/bin/env python3
"""
test_transitions.py

Test the transitions.py module.
"""

import numpy as np
import pytest
from pkg_resources import resource_filename
import amespahdbpythonsuite

from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite import transitions
from amespahdbpythonsuite import data


@pytest.fixture(scope="module")
def pahdb_theoretical():
    xml = 'resources/pahdb-theoretical_cutdown.xml'
    pahdb = AmesPAHdb(filename=resource_filename('amespahdbpythonsuite', xml),
                      check=False, cache=False, update=False)
    return pahdb


@pytest.fixture(scope="module")
def test_spec():

    file1 = 'resources/uid_18_gaussian_6eV_cascade_convolved_test_spec.npy'
    file2 = 'resources/uid_18_drude_6eV_cascade_convolved_test_spec.npy'
    file3 = 'resources/uid_18_lorentzian_6eV_cascade_convolved_test_spec.npy'

    spec1 = np.load(resource_filename('amespahdbpythonsuite', file1))
    spec2 = np.load(resource_filename('amespahdbpythonsuite', file2))
    spec3 = np.load(resource_filename('amespahdbpythonsuite', file3))

    return spec1, spec2, spec3


class TestTransitions():
    """
    Test Transitions class.

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
                                        data=pahdb._AmesPAHdb__getkeybyuids(
                                            'transitions', uids),
                                        pahdb=pahdb._AmesPAHdb__data,
                                        model={'type': 'zerokelvin_m',
                                               'temperature': 0.0,
                                               'description': ''},
                                        units={'abscissa': {'unit': 1,
                                                            'str': 'frequency [wavenumber]'},
                                               'ordinate': {'unit': 2,
                                                            'str': 'integrated cross-section' + '[km/mol]'}})

        test_dict = trans.get()
        key_list = ['type', 'database', 'version',
                    'data', 'uids', 'model', 'units', 'shift']
        assert list(test_dict.keys()) == key_list
        assert test_dict['type'] == 'Transitions'
        assert test_dict['database'] == 'theoretical'
        assert test_dict['version'] == pahdb._AmesPAHdb__data['version']
        # Test set method of data module with provided dictionary.
        test_dict['type'] = 'Data'
        db = data.Data(test_dict, pahdb=pahdb._AmesPAHdb__data)
        assert db.type == 'theoretical'

    def test_getset(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18, 73]
        trans1 = pahdb.gettransitionsbyuid(uids)
        d1 = trans1.get()
        assert(d1['type'] == 'Transitions')
        trans2 = transitions.Transitions()
        trans2.set(d1)
        d2 = trans2.get()
        assert(d2['type'] == 'Transitions')

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
        np.testing.assert_allclose(
            dtest['intensity'], 6.420001406551514e-14)

    def test_calculatedtemperature(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18]
        trans = pahdb.gettransitionsbyuid(uids)
        # Apply calculatedtemperature emission model at 6 eV.
        trans.calculatedtemperature(6 * 1.603e-12)
        # Assert attained Tmax.
        dtest = next(
            (sub for sub in trans.model['temperatures'] if sub['uid'] == uids[0]), None)
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
        dtest_a = next(
            (sub for sub in trans.model['temperatures'] if sub['uid'] == uids[0]), None)
        # Retrieve dictionary at 3068.821 frequency.
        dtest_b = [x for x in trans.data[18] if x['frequency'] == 3068.821][0]
        # Assert attained intensity and temperature.
        assert dtest_a['temperature'] == 1279.7835033561428
        np.testing.assert_allclose(
            dtest_b['intensity'], 1.6710637100014386e-12)

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
        np.testing.assert_allclose(
            test_i['intensity'], 1.6710637100014386e-12)

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
        spec1 = trans.convolve(grid=waven, fwhm=15.0,
                               gaussian=True, multiprocessing=False)
        spec2 = trans.convolve(grid=waven, fwhm=15.0,
                               drude=True, multiprocessing=False)
        spec3 = trans.convolve(grid=waven, fwhm=15.0, multiprocessing=False)
        # Expected gaussian convolved transitions.
        assert spec1.model['temperatures'][0]['temperature'] == 1279.7835033561428
        np.testing.assert_allclose(test_spec[0], spec1.data[18])
        np.testing.assert_allclose(test_spec[1], spec2.data[18])
        np.testing.assert_allclose(test_spec[2], spec3.data[18])

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
                                            width=0.5 * 15.0 /
                                            np.sqrt(2.0 * np.log(2.0)),
                                            x=np.asarray(waven),
                                            gaussian=True,
                                            drude=False,
                                            uid=uids[0])
        np.testing.assert_allclose(test_spec[0], spec[18])

    def test_plot_transitions(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18]
        trans = pahdb.gettransitionsbyuid(uids)
        # Call plot method.
        trans.plot()
