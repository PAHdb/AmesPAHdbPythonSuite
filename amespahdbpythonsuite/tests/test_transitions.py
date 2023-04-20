#!/usr/bin/env python3
"""
test_transitions.py

Test the transitions.py module.
"""

from os.path import exists

import matplotlib.pyplot as plt
import numpy as np
import pytest
from pkg_resources import resource_filename

from amespahdbpythonsuite import transitions
from amespahdbpythonsuite.amespahdb import AmesPAHdb


@pytest.fixture(scope="module")
def pahdb_theoretical():
    xml = "resources/pahdb-theoretical_cutdown.xml"
    db = AmesPAHdb(
        filename=resource_filename("amespahdbpythonsuite", xml),
        check=False,
        cache=False,
        update=False,
    )
    return db


@pytest.fixture(scope="module")
def test_transitions(pahdb_theoretical):
    trans = pahdb_theoretical.gettransitionsbyuid([18])
    trans.cascade(6.0 * 1.603e-12, multiprocessing=False)
    trans.shift(-15.0)
    return trans


@pytest.fixture(scope="module")
def test_spec():
    file1 = "resources/uid_18_gaussian_6eV_cascade_convolved_test_spec.npy"
    file2 = "resources/uid_18_drude_6eV_cascade_convolved_test_spec.npy"
    file3 = "resources/uid_18_lorentzian_6eV_cascade_convolved_test_spec.npy"

    spec1 = np.load(resource_filename("amespahdbpythonsuite", file1))
    spec2 = np.load(resource_filename("amespahdbpythonsuite", file2))
    spec3 = np.load(resource_filename("amespahdbpythonsuite", file3))

    return spec1, spec2, spec3


@pytest.fixture(scope="module")
def test_path(tmp_path_factory):
    d = tmp_path_factory.mktemp("test_transitions")
    return f"{d}/result"


class TestTransitions:
    """
    Test Transitions class.

    """

    def test_instance(self):
        assert isinstance(transitions.Transitions(), transitions.Transitions)

    def test_getset(self, test_transitions):
        t1 = test_transitions.get()
        assert t1["type"] == "Transitions"
        trans2 = transitions.Transitions()
        trans2.set(t1)
        t2 = trans2.get()
        assert t2["type"] == "Transitions"

    def test_set_database_mismatch(self, capsys, pahdb_theoretical):
        transitions.Transitions(
            database="invalid",
            version=pahdb_theoretical.getversion(),
            pahdb=pahdb_theoretical.getdatabaseref(),
        )
        captured = capsys.readouterr()
        assert "DATABASE MISMATCH" in captured.out

    def test_set_version_mismatch(self, capsys, pahdb_theoretical):
        transitions.Transitions(
            database="theoretical",
            version="9999",
            pahdb=pahdb_theoretical.getdatabaseref(),
        )
        captured = capsys.readouterr()
        assert "VERSION MISMATCH" in captured.out

    def test_intersect(self, pahdb_theoretical):
        trans = pahdb_theoretical.gettransitionsbyuid([18, 73, 726, 2054, 223])
        sub_uids = [18, 223]
        trans.intersect(sub_uids)
        assert list(trans.uids) == sub_uids
        assert list(trans.data.keys()) == sub_uids

    def test_difference(self, pahdb_theoretical):
        sub_uids = [18, 223]
        trans = pahdb_theoretical.gettransitionsbyuid([18, 73, 726, 2054, 223])
        trans.difference(sub_uids)
        assert list(trans.uids) == [73, 726, 2054]
        assert list(trans.data.keys()) == [73, 726, 2054]

    def test_shift(self, pahdb_theoretical):
        trans = pahdb_theoretical.gettransitionsbyuid([18])
        dtest = [x for x in trans.data[18] if x["frequency"] == 3068.821][0]
        trans.shift(-15.0)
        assert dtest["frequency"] == 3053.821

    def test_fixedtemperature(self, pahdb_theoretical):
        trans = pahdb_theoretical.gettransitionsbyuid([18])
        trans.fixedtemperature(600)
        dtest = [x for x in trans.data[18] if x["frequency"] == 3068.821][0]
        np.testing.assert_allclose(dtest["intensity"], 6.420001406551514e-14)

    def test_calculatedtemperature(self, pahdb_theoretical):
        trans = pahdb_theoretical.gettransitionsbyuid([18])
        trans.calculatedtemperature(6 * 1.603e-12)
        dtest = next(
            (sub for sub in trans.model["temperatures"] if sub["uid"] == 18), None
        )
        assert dtest["temperature"] == 1279.7835033561428

    def test_cascade(self, pahdb_theoretical):
        trans = pahdb_theoretical.gettransitionsbyuid([18])
        trans.cascade(6 * 1.603e-12, multiprocessing=False)
        dtest_a = next(
            (sub for sub in trans.model["temperatures"] if sub["uid"] == 18), None
        )
        dtest_b = [x for x in trans.data[18] if x["frequency"] == 3068.821][0]
        assert dtest_a["temperature"] == 1279.7835033561428
        np.testing.assert_allclose(dtest_b["intensity"], 1.6710637100014386e-12)

    def test_partial_cascade(self, pahdb_theoretical):
        trans_multi = pahdb_theoretical.gettransitionsbyuid([18])
        data = trans_multi.get()
        intf, tmax = transitions.Transitions._cascade_em_model(6 * 1.603e-12, data['data'][18])
        test_i = [x for x in intf if x["frequency"] == 3068.821][0]
        assert tmax == 1279.7835033561428
        np.testing.assert_allclose(test_i["intensity"], 1.6710637100014386e-12)

    def test_convolve_gaussian(self, test_transitions, test_spec):
        spec = test_transitions.convolve(
            grid=[1e4 / x for x in np.arange(5, 20, 0.4)],
            fwhm=15.0,
            gaussian=True,
            multiprocessing=False,
        )
        np.testing.assert_allclose(test_spec[0], spec.data[18])

    def test_convolve_drude(self, test_transitions, test_spec):
        spec = test_transitions.convolve(
            grid=[1e4 / x for x in np.arange(5, 20, 0.4)],
            fwhm=15.0,
            drude=True,
            multiprocessing=False,
        )
        np.testing.assert_allclose(test_spec[1], spec.data[18])

    def test_convolve_lorentzian(self, test_transitions, test_spec):
        spec = test_transitions.convolve(
            grid=[1e4 / x for x in np.arange(5, 20, 0.4)],
            fwhm=15.0,
            multiprocessing=False,
        )
        np.testing.assert_allclose(test_spec[2], spec.data[18])

    def test_partial_convolve(self, pahdb_theoretical, test_spec):
        trans_multi = pahdb_theoretical.gettransitionsbyuid([18])
        k = [1e4 / x for x in np.arange(5, 20, 0.4)]
        trans_multi.cascade(6.0 * 1.603e-12, multiprocessing=False)
        trans_multi.shift(-15.0)
        data = trans_multi.get()
        spec = transitions.Transitions._get_intensities(
            npoints=len(k),
            xmin=np.min(k),
            xmax=np.max(k),
            clip=3,
            width=0.5 * 15.0 / np.sqrt(2.0 * np.log(2.0)),
            x=np.asarray(k),
            gaussian=True,
            drude=False,
            data=data["data"][18],
        )
        np.testing.assert_allclose(test_spec[0], spec)

    def test_plot_transitions(self, monkeypatch, test_transitions):
        monkeypatch.setattr(plt, "show", lambda: None)
        test_transitions.plot(Show=True)

    def test_write_transitions(self, test_transitions, test_path):
        test_transitions.write(f"{test_path}.tbl")
        assert exists(f"{test_path}.tbl")
