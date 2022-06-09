#!/usr/bin/env python3
"""
test_species.py

Test the species.py module.
"""

import pytest
from pkg_resources import resource_filename

from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite import geometry, species
from amespahdbpythonsuite import transitions
from amespahdbpythonsuite import laboratory


@pytest.fixture(scope="module")
def species_test():
    xml = 'resources/pahdb-theoretical_cutdown.xml'
    db = AmesPAHdb(filename=resource_filename('amespahdbpythonsuite', xml),
                   check=False, cache=False, update=False)
    s = db.getspeciesbyuid([18, 73])
    return s


class TestSpecies():
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
        assert(s1['type'] == 'Species')
        species2 = species.Species()
        species2.set(s1)
        s2 = species2.get()
        assert(s2['type'] == 'Species')
