#!/usr/bin/env python3
"""
test_amespahdb.py

Test the amespahdb.py module.
"""

from pkg_resources import resource_filename
import amespahdbpythonsuite

from amespahdbpythonsuite.amespahdb import AmesPAHdb


class TestAmesPAHdb():
    """
    Test AmesPAHdb class.

    """
    def test_instance(self):
        # Define database name and location.
        xml = 'resources/pahdb-theoretical_cutdown.xml'
        database = resource_filename('amespahdbpythonsuite', xml)
        # Read the database.
        pahdb = AmesPAHdb(filename=database, check=False, cache=True)
        assert isinstance(pahdb, amespahdbpythonsuite.amespahdb.AmesPAHdb)

    def test_keybyuids(self):
        # Define database name and location.
        xml = 'resources/pahdb-theoretical_cutdown.xml'
        database = resource_filename('amespahdbpythonsuite', xml)
        # Read the database.
        pahdb = AmesPAHdb(filename=database, check=False, cache=True)
        # UIDs test list.
        uids = [18, 73, 726, 2054, 223]
        res = pahdb._AmesPAHdb__getkeybyuids('formula', uids)
        assert sorted(uids) == sorted(list(res.keys()))

    def test_gettransitionsbyuid(self):
        # Define database name and location.
        xml = 'resources/pahdb-theoretical_cutdown.xml'
        database = resource_filename('amespahdbpythonsuite', xml)
        # Read the database.
        pahdb = AmesPAHdb(filename=database, check=False, cache=True)
        # UIDs test list.
        uids = [18, 73, 726, 2054, 223]
        trans = pahdb.gettransitionsbyuid(uids)
        assert isinstance(trans, amespahdbpythonsuite.transitions.Transitions)

    def test_getlaboratorybyuid(self):
        # UIDs test list.
        uids = [273]
        # Retrieve database dictionary of provided UIDs.
        xml = 'resources/pahdb-experimental_cutdown.xml'
        database = resource_filename('amespahdbpythonsuite', xml)
        pahdb = AmesPAHdb(filename=database, check=False, cache=True)
        lab = pahdb.getlaboratorybyuid(uids)
        assert isinstance(lab, amespahdbpythonsuite.laboratory.Laboratory)

    def test_getspeciesbyuid(self):
        # Define database name and location.
        xml = 'resources/pahdb-theoretical_cutdown.xml'
        database = resource_filename('amespahdbpythonsuite', xml)
        # Read the database.
        pahdb = AmesPAHdb(filename=database, check=False, cache=True)
        # UIDs test list.
        uids = [18, 726, 2054]
        spec = pahdb.getspeciesbyuid(uids)
        assert isinstance(spec, amespahdbpythonsuite.species.Species)

    def test_getgeometrybyuid(self):
        # Define database name and location.
        xml = 'resources/pahdb-theoretical_cutdown.xml'
        database = resource_filename('amespahdbpythonsuite', xml)
        # Read the database.
        pahdb = AmesPAHdb(filename=database, check=False, cache=True)
        # UIDs test list.
        uids = [18, 73, 726, 2054, 223]
        geo = pahdb.getgeometrybyuid(uids)
        assert isinstance(geo, amespahdbpythonsuite.geometry.Geometry)
