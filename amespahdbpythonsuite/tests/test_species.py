#!/usr/bin/env python3
"""
test_species.py

Test the species.py module.
"""

import amespahdbpythonsuite

from amespahdbpythonsuite.species import Species


class TestSpecies():
    """
    Test Species class.

    """
    def test_instance(self):
        spec = Species()
        assert isinstance(spec, amespahdbpythonsuite.species.Species)
