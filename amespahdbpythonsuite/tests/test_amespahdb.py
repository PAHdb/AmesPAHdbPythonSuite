#!/usr/bin/env python3
"""
test_amespahdb.py

Test the amespahdb.py module.
"""
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


@pytest.fixture(scope="module")
def pahdb_laboratory():
    xml = 'resources/pahdb-experimental_cutdown.xml'
    pahdb = AmesPAHdb(filename=resource_filename('amespahdbpythonsuite', xml),
                      check=False, cache=False)
    return pahdb


class TestAmesPAHdb():
    """
    Test AmesPAHdb class.

    """

    def test_instance(self, pahdb_theoretical):
        # Read the database.
        assert isinstance(pahdb_theoretical,
                          amespahdbpythonsuite.amespahdb.AmesPAHdb)

    def test_env(self):
        # TODO: Turn the sys.exit into exceptions.
        import os
        if 'AMESPAHDEFAULTDB' in os.environ:
            del os.environ["AMESPAHDEFAULTDB"]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            AmesPAHdb()
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 1

    def test_file_not_exist(self):
        # TODO: Turn the sys.exit into exceptions.
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            AmesPAHdb(filename='file_does_not_exist.xml')
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 2

    def test_cache(self, capsys):
        xml = 'resources/pahdb-experimental_cutdown.xml'
        # Make sure the AmesPAHdb module is called to create the cached database.
        with capsys.disabled():
            AmesPAHdb(filename=resource_filename('amespahdbpythonsuite', xml),
                      check=False, cache=True)
        AmesPAHdb(filename=resource_filename('amespahdbpythonsuite', xml),
                  check=False, cache=True)
        capture = capsys.readouterr()
        assert capture.out.find('RESTORING DATABASE FROM CACHE') >= 0

    def test_checkversion(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        assert not pahdb.checkversion('wrong')

    def test_keybyuids(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18, 73, 726, 2054, 223]
        res = pahdb._AmesPAHdb__getkeybyuids('formula', uids)
        assert sorted(uids) == sorted(list(res.keys()))

    def test_gettransitionsbyuid(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18, 73, 726, 2054, 223]
        trans = pahdb.gettransitionsbyuid(uids)
        assert isinstance(trans, amespahdbpythonsuite.transitions.Transitions)

    def test_getlaboratorybyuid(self, pahdb_laboratory):
        # Read the database.
        pahdb = pahdb_laboratory
        # UIDs test list.
        uids = [273]
        lab = pahdb.getlaboratorybyuid(uids)
        assert isinstance(lab, amespahdbpythonsuite.laboratory.Laboratory)

    def test_getspeciesbyuid(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18, 726, 2054]
        spec = pahdb.getspeciesbyuid(uids)
        assert isinstance(spec, amespahdbpythonsuite.species.Species)

    def test_getgeometrybyuid(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18, 73, 726, 2054, 223]
        geo = pahdb.getgeometrybyuid(uids)
        assert isinstance(geo, amespahdbpythonsuite.geometry.Geometry)

    def test_getdatabaseref(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        assert isinstance(pahdb.getdatabaseref(), dict)

    def test_search_formula(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        uids = pahdb.search("c24h12")
        assert uids == [18]

    def test_search_charge(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        uids = pahdb.search("charge=0")
        assert uids == [18, 73, 726, 2054, 223]

    def test_search_size_and_charge(self, pahdb_theoretical):
        # Read the database.
        pahdb = pahdb_theoretical
        uids = pahdb.search("c=24 and h=12 and neutral")
        assert uids == [18]
