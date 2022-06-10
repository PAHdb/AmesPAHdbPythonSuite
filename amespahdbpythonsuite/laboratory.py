#!/usr/bin/env python3

from amespahdbpythonsuite.data import Data


class Laboratory(Data):
    """
    AmesPAHdbPythonSuite laboratory class.
    Contains methods to work with a laboratory spectrum.

    """

    def __init__(self, d=None, **keywords):
        super().__init__(d, **keywords)
        self.set(d, **keywords)

    def set(self, d, **keywords):
        """
        Calls class: :class:`amespahdbpythonsuite.spectrum.Spectrum.set to parse keywords.

        """
        Data.set(self, d, **keywords)

    def get(self):
        """
        Calls class: :class:`amespahdbpythonsuite.spectrum.Spectrum.get.
        Assigns class variables from inherited dictionary.

        """
        d = Data.get(self)
        d['type'] = self.__class__.__name__

        return d
