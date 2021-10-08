#!/usr/bin/env python3
"""
test_amespahdb.py

Test the Transitions.py module.
"""

import amespahdbpythonsuite

from amespahdbpythonsuite.transitions import Transitions


class TestTransitions():
    """
    Test AmesPAHdb class.

    """
    def test_instance(self):
        trans = Transitions()
        assert isinstance(trans, amespahdbpythonsuite.transitions.Transitions)
