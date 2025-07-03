#!/usr/bin/env python3
"""
test_species.py

Test the species.py module.
"""

import copy

import importlib_resources
import pytest

from amespahdbpythonsuite import geometry, laboratory, species, transitions
from amespahdbpythonsuite.amespahdb import AmesPAHdb


@pytest.fixture(scope="module")
def species_test():
    xml = (
        importlib_resources.files("amespahdbpythonsuite")
        / "resources/pahdb-theoretical_cutdown.xml"
    )
    db = AmesPAHdb(
        filename=xml,
        check=False,
        cache=False,
        update=False,
    )
    s = db.getspeciesbyuid([18, 73, 726, 2054, 223])
    return s


class TestSpecies:
    """
    Test Species class.

    """

    def test_instance(self):
        assert isinstance(species.Species(), species.Species)

    def test_transitions(self, species_test):
        assert isinstance(species_test.transitions(), transitions.Transitions)

    def geometry(self, species_test):
        assert isinstance(species_test.geometry(), geometry.Geometry)

    def laboratory(self, species_test):
        assert isinstance(species_test.laboratory(), laboratory.Laboratory)

    def test_references(self, species_test):
        assert isinstance(species_test.references(), dict)

    def test_comments(self, species_test):
        assert isinstance(species_test.comments(), dict)

    def test_getset(self, species_test):
        s1 = species_test.get()
        assert s1["type"] == "Species"
        species2 = species.Species()
        species2.set(s1)
        s2 = species2.get()
        assert s2["type"] == "Species"

    def test_getset_mismatch_database(self, capsys, species_test):
        s = copy.copy(species_test)
        s.set(database="none")
        capture = capsys.readouterr()
        assert capture.out.find("DATABASE MISMATCH")

    def test_getset_mismatch_version(self, capsys, species_test):
        s = copy.copy(species_test)
        s.set(version="none")
        capture = capsys.readouterr()
        assert capture.out.find("VERSION MISMATCH")

    def test_intersect(self, species_test):
        sub_uids = [18, 223]
        s = copy.copy(species_test)
        s.intersect(sub_uids)
        assert list(s.uids) == sub_uids
        assert list(s.data.keys()) == sub_uids

    def test_difference(self, species_test):
        sub_uids = [18, 223]
        s = copy.copy(species_test)
        s.difference(sub_uids)
        assert list(s.uids) == [73, 726, 2054]
        assert list(s.data.keys()) == [73, 726, 2054]

    def test_formatformula(self):
        assert (
            species.formatformula("C10H8++")
            == r"C$_{\mathregular{10}}H$_{\mathregular{8}}$^{\mathregular{++}}"
        )

    def test_print(self, species_test):
        assert species_test.print(18, str=True)[0:14] == "UID       : 18"
        assert species_test.print(str=True)[0:14] == "=============="
