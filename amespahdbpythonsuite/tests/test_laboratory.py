#!/usr/bin/env python3
"""
test_amespahdb.py

Test the Transitions.py module.
"""

import amespahdbpythonsuite

from amespahdbpythonsuite.laboratory import Laboratory


class TestTransitions():
    """
    Test AmesPAHdb class.

    """
    def test_instance(self):
        lab = Laboratory()
        assert isinstance(lab, amespahdbpythonsuite.laboratory.Laboratory)
