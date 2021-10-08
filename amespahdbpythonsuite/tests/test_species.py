#!/usr/bin/env python3
"""
test_amespahdb.py

Test the Transitions.py module.
"""

import amespahdbpythonsuite

from amespahdbpythonsuite.species import Species


class TestTransitions():
    """
    Test AmesPAHdb class.

    """
    def test_instance(self):
        spec = Species()
        assert isinstance(spec, amespahdbpythonsuite.species.Species)
