#!/usr/bin/env python3

from amespahdbpythonsuite.spectrum import Spectrum


class Fitted(Spectrum):
    """
    AmesPAHdbPythonSuite fitted class

    """

    def __init__(self, d=None, **keywords):
        Spectrum.__init__(self, d, **keywords)
