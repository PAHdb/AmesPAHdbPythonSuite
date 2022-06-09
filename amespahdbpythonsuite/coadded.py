#!/usr/bin/env python3

from amespahdbpythonsuite.spectrum import Spectrum


class Coadded(Spectrum):
    """
    AmesPAHdbPythonSuite coadded class

    """

    def __init__(self, d=None, **keywords):
        Spectrum.__init__(self, d, **keywords)
        self.weights = []
        self.averaged = False
        self.__set(d, **keywords)
        return None

    def set(self, d, **keywords):
        """
        Calls class: :class:`amespahdbpythonsuite.Spectrum.spectrum.set to parse keywords.

        """
        Spectrum.set(self, d, **keywords)
        self.__set(d, **keywords)

    def __set(self, d, **keywords):
        """
        Populate data dictionary helper.

        """
        if d:
            if d.get('type', '') == self.__class__.__name__:
                if 'weights' not in keywords:
                    self.weights = d['weights']
                if 'averaged' not in keywords:
                    self.averaged = d['averaged']

        if len(keywords.get('weights', [])):
            self.weights = keywords.get('weights')
        if 'averaged' in keywords:
            self.averaged = keywords.get('averaged')

    def get(self):
        """
        Calls class: :class:`amespahdbpythonsuite.spectrum.Spectrum.get.
        Assigns class variables from inherited dictionary.

        """
        d = Spectrum.get(self)
        d['type'] = self.__class__.__name__
        d['weights'] = self.weights
        d['averaged'] = self.averaged

        return d
