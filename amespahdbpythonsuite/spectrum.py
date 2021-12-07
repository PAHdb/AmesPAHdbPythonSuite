#!/usr/bin/env python3

import numpy as np
import copy

from scipy import optimize

from amespahdbpythonsuite.transitions import Transitions


class Spectrum(Transitions):
    """
    AmesPAHdbPythonSuite spectrum class.
    Contains methods to fit and plot the input spectrum.

    """

    def __init__(self, d=None, **keywords):
        Transitions.__init__(self, d, **keywords)

    def set(self, d=None, **keywords):
        """
        Calls class: :class:`amespahdbpythonsuite.transitions.Transitions.set to parse keywords.

        """
        Transitions.set(self, d, **keywords)
        if d:
            if d.get('type', '') == self.__class__.__name__:
                if not keywords.get('grid'):
                    self.grid = d['grid']
                if not keywords.get('profile'):
                    self.profile = d['profile']
                if not keywords.get('fwhm'):
                    self.fwhm = d['fwhm']

                Transitions.set(self, d, **keywords)

        if len(keywords.get('grid', [])):
            self.grid = keywords.get('grid')
        if keywords.get('profile'):
            self.profile = keywords.get('profile')
        if keywords.get('fwhm'):
            self.fwhm = keywords.get('fwhm')

    def get(self):
        """
        Calls class: :class:`amespahdbpythonsuite.transitions.Transitions.get.
        Assigns class variables from inherited dictionary.

        """
        d = Transitions.get(self)
        d['type'] = self.__class__.__name__
        d['grid'] = self.grid
        d['profile'] = self.profile
        d['fwhm'] = self.fwhm

        return d

    def subdatabase(self, uids):
        """
        Creates a subset of the database containing the retrieved UIDs.
        """
        d = {}
        for key in self.pahdb:
            if key == 'species':
                d[key] = dict((uid, self.pahdb['species'][uid]) for uid in uids)
            else:
                d[key] = self.pahdb[key]
        return copy.deepcopy(d)

    def fit(self, y, yerr=None, **keywords):
        """
        Fits the input spectrum.

        """

        matrix = []
        matrix = np.array(list(self.data.values()))

        if yerr is None:
            # Do NNLS.
            b = list(y)
            m = matrix
        else:
            # Do NNLC.
            b = list(np.divide(y, yerr))
            m = np.divide(matrix, yerr)

        solution, norm = optimize.nnls(m.T, b)

        # Initialize lists and dictionaries.
        uids = []
        data = {}
        weights = {}

        # Retrieve uids, data, and fit weights dictionaries.
        for uid, s, m, in zip(self.uids, solution, matrix):
            if s > 0:
                intensities = []
                uids.append(uid)
                for d in m:
                    intensities.append(s * d)
                data[uid] = np.array(intensities)
                weights[uid] = s

        # Call Fitted Class to plot the fitted spectrum and get the fit breakdown.
        from amespahdbpythonsuite.fitted import Fitted

        return Fitted(type=self.type,
                      version=self.version,
                      data=data,
                      pahdb=self.subdatabase(uids),
                      uids=uids,
                      model=self.model,
                      units=self.units,
                      shift=self.shift,
                      grid=self.grid,
                      profile=self.profile,
                      fwhm=self.fwhm,
                      observation=list(y),
                      weights=weights)

    def plot(self, **keywords):
        """
        Plot the spectrum.

        """
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm

        _, ax = plt.subplots()
        ax.minorticks_on()
        ax.tick_params(which='major', right='on', top='on', direction='in', length=5)
        ax.tick_params(which='minor', right='on', top='on', direction='in', length=3)
        colors = cm.rainbow(np.linspace(0, 1, len(self.uids)))

        for uid, col in zip(self.uids, colors):
            y = self.data[uid]
            ax.plot(self.grid, y, color=col)

        ax.set_xlim((max(self.grid), min(self.grid)))
        ax.tick_params(axis='both', labelsize=14)

        ax.set_xlabel(self.units['abscissa']['str'], fontsize=14)
        ax.set_ylabel(self.units['ordinate']['str'], fontsize=14)

        plt.show()
