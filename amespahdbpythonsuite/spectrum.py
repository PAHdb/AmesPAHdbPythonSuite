#!/usr/bin/env python3

from amespahdbpythonsuite.transitions import Transitions


class Spectrum(Transitions):
    """
    AmesPAHdbPythonSuite spectrum class

    """

    def __init__(self, d=None, **keywords):
        Transitions.__init__(self, d, **keywords)
        return None
