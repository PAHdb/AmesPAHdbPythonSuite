#!/usr/bin/env python3

import operator
import numpy as np

from scipy import integrate
from astropy.io import ascii

from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite.spectrum import Spectrum

message = AmesPAHdb.message


class Fitted(Spectrum):
    """
    AmesPAHdbPythonSuite fitted class

    """

    def __init__(self, d=None, **keywords):
        """
        Initialize fitted class.

        """
        Spectrum.__init__(self, d, **keywords)
        self.weights = None
        self.observation = None
        self.atoms = None
        self.__set(d, **keywords)

    def plot(self, **keywords):
        """
        Plotting method for the fitted spectrum and breakdown components.

        """
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gs
        import matplotlib.cm as cm

        mpl.rc('xtick', labelsize=12)
        mpl.rc('ytick', labelsize=12)

        if keywords.get('residual', False):
            figures = plt.figure()
            spec = gs.GridSpec(2, 1, height_ratios=[3, 1])
            axis = plt.subplot(spec[0])
            axis = [axis, plt.subplot(spec[1], sharex=axis)]
            figures.subplots_adjust(hspace=0)
            axis[0].tick_params(axis='x', which='both', labelbottom='off')
        elif keywords.get('charge', False) \
                or keywords.get('size', False) \
                or keywords.get('composition', False):
            figures, axis = plt.subplots()
            axis = [axis]
        else:
            figures = plt.figure()
            spec = gs.GridSpec(1, 2, width_ratios=[2, 3])
            axis = plt.subplot(spec[0])
            axis = [axis, plt.subplot(spec[1])]
            figures.subplots_adjust(wspace=0.25)
            axis[0].tick_params(axis='x', which='both',
                                bottom='off', top='off', labelbottom='off')
            axis[0].tick_params(axis='y', which='both',
                                left='off', right='off', labelleft='off')
            axis[0].set_ylim((0, 1))
            axis = list(reversed(axis))

        if keywords.get('wavelength', False):
            x = 1e4 / self.grid
            axis[0].set_xlim((min(x), max(x)))
            xtitle = 'Wavelength'
        else:
            x = self.grid
            axis[0].set_xlim((max(x), min(x)))
            xtitle = self.units['abscissa']['str']

        axis[0].errorbar(x, self.observation,
                         yerr=keywords['sigma'],
                         fmt='o', mfc='white', color='k', ecolor='k',
                         markersize=3, elinewidth=0.2, capsize=0.8, label='obs')
        if keywords.get('title', False):
            axis[0].set_title(keywords['title'])

        axis[0].minorticks_on()
        axis[0].tick_params(which='major', right='on',
                            top='on', direction='in', length=5)
        axis[0].tick_params(which='minor', right='on',
                            top='on', direction='in', length=3)

        colors = cm.rainbow(np.linspace(0, 1, len(self.uids)))
        for uid, col in zip(self.uids, colors):
            y = self.data[uid]
            if not keywords.get('charge', False) \
                    and not keywords.get('size', False) \
                    and not keywords.get('composition', False):
                axis[0].plot(x, y, color=col)

        axis[0].plot(x, self.getfit(), color='tab:purple', label='fit')

        if keywords.get('charge', False) \
                or keywords.get('size', False) \
                or keywords.get('composition', False):
            classes = self.getclasses()
            if keywords.get('charge', False):
                if not isinstance(classes['anion'], int):
                    axis[0].plot(x, classes['anion'], 'r-', label='anion')
                if not isinstance(classes['neutral'], int):
                    axis[0].plot(x, classes['neutral'], 'g-', label='neutral')
                if not isinstance(classes['cation'], int):
                    axis[0].plot(x, classes['cation'], 'b-', label='cation')
                axis[0].axhline(0, linestyle='--', color='gray')
                axis[0].legend(fontsize=12)

            elif keywords.get('size', False):
                if not isinstance(classes['small'], int):
                    axis[0].plot(x, classes['small'], 'r-', label='small')
                if not isinstance(classes['large'], int):
                    axis[0].plot(x, classes['large'], 'g-', label='large')
                axis[0].axhline(0, linestyle='--', color='gray')
                axis[0].legend(fontsize=12)

            elif keywords.get('composition', False):
                if not isinstance(classes['pure'], int):
                    axis[0].plot(x, classes['pure'], 'r-', label='pure')
                if not isinstance(classes['nitrogen'], int):
                    axis[0].plot(x, classes['nitrogen'],
                                 'g-', label='nitrogen')
                axis[0].axhline(0, linestyle='--', color='gray')
                axis[0].legend(fontsize=12)

        elif keywords.get('residual', False):
            y = self.getresidual()
            axis[1].plot(x, y)

        else:
            axis[1].text(0.05, 0.95, ("%s" + 5 * ' ' + "%s")
                         % ('UID', 'WEIGHT'), family='monospace')
            axis[1].xaxis.set_visible(False)
            axis[1].yaxis.set_visible(False)
            ypos = 0.95 - 0.05
            for uid, w, col in zip(self.uids, self.weights.values(), colors):
                str = ("%d" + (5 * ' ') + "%.2e") % (uid, w)
                axis[1].text(0.05, ypos, str, color=col, family='monospace')
                ypos -= 0.05
                if ypos <= 0.05:
                    axis[1].text(0.05, ypos, 'more...', family='monospace')
                    break

        axis[0].set_ylabel(keywords['units'][1], fontsize=14)
        if keywords.get('residual', False):
            axis[1].set_xlabel(
                f"{xtitle} ({keywords['units'][0]})", fontsize=14)
            axis[1].set_ylabel('residual', fontsize=14)
            axis[1].minorticks_on()
            axis[1].tick_params(which='major', right='on',
                                top='on', direction='in', length=5)
            axis[1].tick_params(which='minor', right='on',
                                top='on', direction='in', length=3)

        else:
            axis[0].set_xlabel(
                f"{xtitle} ({keywords['units'][0]})", fontsize=14)
        if keywords.get('show', False):
            plt.show()

        # save plots
        figures.tight_layout()
        if keywords['ftype']:
            figures.savefig(
                f"{keywords['outputname']}_{keywords['ptype']}.{keywords['ftype']}")
        plt.close(figures)

    def set(self, d=None, **keywords):
        """
        Calls class: :class:`amespahdbpythonsuite.spectrum.Spectrum.set to parse keywords.

        """
        Spectrum.set(self, d, **keywords)
        self.__set(d, **keywords)

    def __set(self, d=None, **keywords):
        """
        Populate data dictionary helper.

        """
        if d:
            if d.get('type', '') == self.__class__.__name__:
                if not keywords.get('observation'):
                    self.observation = d['observation']
                if not keywords.get('weights'):
                    self.weights = d['weights']

        if len(keywords.get('observation', [])):
            self.observation = keywords.get('observation')
        if keywords.get('weights'):
            self.weights = keywords.get('weights')

    def get(self):
        """
        Calls class: :class:`amespahdbpythonsuite.spectrum.Spectrum.get.
        Assigns class variables from inherited dictionary.

        """
        d = Spectrum.get(self)
        d['type'] = self.__class__.__name__
        d['observation'] = self.observation
        d['weights'] = self.weights

        return d

    def getchisquared(self):
        """
        Calculates the chi-squared of the fit.

        """
        if self.observation:
            return np.sum(
                (self.observation.data.y - self.getfit()) ** 2 / self.observation.data.ystdev)

        return None

    def getnorm(self):
        """
        Retrieves the Norm of the fit.

        """
        return np.total((self.observation.data.y - self.getfit()) ** 2)

    def getobservation(self):
        """
        Retrieves the ordinate values of the observation.

        """
        return self.observation

    def getfit(self):
        """
        Retrieves the sum of fit values.

        """
        fitsum = sum(self.data.values())
        return fitsum

    def getresidual(self):
        """
        Retrieves the residual of the fit.
        """
        return self.observation - self.getfit()

    def getweights(self):
        """
        Retrieves the weights of the fitted PAHs.

        """
        return self.weights

    def write(self, prefix=''):
        """
        Retrieve fitted PAH properties and write to file.

        """

        # Retrieve properties
        fweight = list(self.weights.values())
        uids = self.uids
        formula = [self.pahdb['species'][uid]['formula'] for uid in self.uids]
        nc = [self.atoms[uid]['nc'] for uid in self.uids]
        charge = [self.pahdb['species'][uid]['charge'] for uid in self.uids]
        mweight = [self.pahdb['species'][uid]['weight'] for uid in self.uids]
        nsolo = [self.pahdb['species'][uid]['n_solo'] for uid in self.uids]
        nduo = [self.pahdb['species'][uid]['n_duo'] for uid in self.uids]
        ntrio = [self.pahdb['species'][uid]['n_trio'] for uid in self.uids]
        nquartet = [self.pahdb['species'][uid]['n_quartet']
                    for uid in self.uids]
        nquintet = [self.pahdb['species'][uid]['n_quintet']
                    for uid in self.uids]

        if prefix == '':
            prefix = self.__class__.__name__

        # Write to file.
        filename = f"{prefix}_results.txt"
        ascii.write([uids, formula, nc, charge, mweight,
                    nsolo, nduo, ntrio, nquartet, nquintet, fweight],
                    filename,
                    names=['#UID', 'formula', 'Nc', 'charge', 'mweight',
                           'n_solo', 'n_duo', 'n_trio', 'n_quartet', 'n_quintet', 'fweight'],
                    overwrite=True)

        message(f'WRITTEN: {filename}')

    def getclasses(self, **keywords):
        """
        Retrieves the spectra of the different classes of the fitted PAHs.

        """
        if not self.pahdb:
            message('VALID DATABASE NEEDED TO GET CLASSES')

            return None

        # Set atom keywords in species dictionary.
        self._atoms()

        # Set subclasses dictionary.
        subclasses = self._subclasses(**keywords)

        classes = dict()

        for key in subclasses:
            classes[key] = self.__classes(subclasses[key])

        puids = [uid for uid in self.uids
                 if self.atoms[uid]['nn'] == 0
                 and self.atoms[uid]['no'] == 0
                 and self.atoms[uid]['nmg'] == 0
                 and self.atoms[uid]['nsi'] == 0
                 and self.atoms[uid]['nfe'] == 0]
        classes['pure'] = sum(
            {k: v for k, v in self.data.items() if k in puids}.values())

        return classes

    def __classes(self, s):
        """
        Retrieves the intensities of a given subclass.

        """

        if s['subclass'] == 'charge':
            uids = [uid for uid in self.uids
                    if s['operator'](self.pahdb['species'][uid][s['subclass']], s['operand'])]
        else:
            uids = [uid for uid in self.uids
                    if s['operator'](self.atoms[uid][s['subclass']], s['operand'])]

        intensities = sum(
            {k: v for k, v in self.data.items() if k in uids}.values())
        return intensities

    def getbreakdown(self, **keywords):
        """
        Retrieves the breakdown of the fitted PAHs.

        """
        if not self.pahdb:
            message('VALID DATABASE NEEDED TO GET CLASSES')

            return None

        # Set atom dictionary if it doesn't exist.
        if not self.atoms:
            self._atoms()

        # Getting fit weights
        fweight = np.array(list(self.weights.values()))

        breakdown = {'solo': np.sum([self.pahdb['species'][uid]['n_solo'] for uid in self.uids]),
                     'duo': np.sum([self.pahdb['species'][uid]['n_duo'] for uid in self.uids]),
                     'trio': np.sum([self.pahdb['species'][uid]['n_trio'] for uid in self.uids]),
                     'quartet': np.sum([self.pahdb['species'][uid]['n_quartet'] for uid in self.uids]),
                     'quintet': np.sum([self.pahdb['species'][uid]['n_quintet'] for uid in self.uids])
                     }

        total = 1.0

        if keywords.get('flux', False):
            if not self.grid:
                message('GRID IS NOT SET')

                return None

            classes = self.getclasses(**keywords)

            if not keywords.get('absolute', False):
                total, err = integrate.quad(self.grid, self.getfit())
                print('first integrate occurance')

            for key in classes:
                breakdown[key] = integrate.quad(
                    self.grid, classes[key]) / total
                print('second integrate occurance')

            return breakdown

        if not self.weights:
            message('WEIGHTS ARE NOT SET')

            return None

        if not keywords.get('absolute', False):
            total = np.sum(fweight)

        # Set subclasses dictionary.
        subclasses = self._subclasses(**keywords)

        for key in subclasses:
            breakdown[key] = self.__breakdown(subclasses[key]) / total

        # Obtaining pure PAH breakdown.
        uids = [uid for uid in self.uids
                if self.atoms[uid]['nn'] == 0
                and self.atoms[uid]['no'] == 0
                and self.atoms[uid]['nmg'] == 0
                and self.atoms[uid]['nsi'] == 0
                and self.atoms[uid]['nfe'] == 0]

        if len(uids) > 0:
            breakdown['pure'] = np.sum([self.weights[uid]
                                       for uid in uids]) / total

        # Getting Nc.
        nc = np.array([self.atoms[uid]['nc'] for uid in self.uids])
        breakdown['nc'] = np.sum(nc * fweight) / np.sum(fweight)
        # Getting fit uncertainty.
        pahdberr = self.__geterror()
        for key in pahdberr.keys():
            breakdown[key] = pahdberr[key]

        return breakdown

    def __breakdown(self, s):
        """
        Retrieve the sum of the fitting weights for the fitted PAHs.

        """
        if s['subclass'] == 'charge':
            uids = [uid for uid in self.uids
                    if s['operator'](self.pahdb['species'][uid][s['subclass']], s['operand'])]
        else:
            uids = [uid for uid in self.uids
                    if s['operator'](self.atoms[uid][s['subclass']], s['operand'])]

        if len(uids) > 0:
            return np.sum([self.weights[uid] for uid in uids])

        else:
            return 0.0

    def _subclasses(self, **keywords):
        """
        Create subclasses dictionary.

        """

        subclasses = {'anion': {'subclass': 'charge',
                                'operator': operator.lt,
                                'operand': 0},
                      'neutral': {'subclass': 'charge',
                                  'operator': operator.eq,
                                  'operand': 0},
                      'cation': {'subclass': 'charge',
                                 'operator': operator.gt,
                                 'operand': 0},
                      'small': {'subclass': 'nc',
                                'operator': operator.le,
                                'operand': keywords.get('small', 50)},
                      'large': {'subclass': 'nc',
                                'operator': operator.gt,
                                'operand': keywords.get('small', 50)},
                      'nitrogen': {'subclass': 'nn',
                                   'operator': operator.gt,
                                   'operand': 0}}

        return subclasses

    def _atoms(self):
        """
        Create atoms dictionary.

        """

        # Create reference dictionary with atomic numbers for c, h, n, o, mg, si, and fe.
        nelem = {'nc': 6, 'nh': 1, 'nn': 7, 'no': 8,
                 'nmg': 12, 'nsi': 14, 'nfe': 26}

        # Initialize atoms dictionary.
        self.atoms = {key: {} for key in self.uids}

        for uid in self.uids:
            # Initialize dictionary based on reference dictionary.
            dnelem = dict.fromkeys(nelem)
            # Loop through the keys to determine the number of each atom present in a given uid.
            for key in dnelem.keys():
                dnelem[key] = len([sub['type'] for sub in self.pahdb['species'][uid]['geometry']
                                  if sub['type'] == nelem[key]])
            self.atoms[uid] = dnelem

    def __geterror(self):
        """
        Calculates the PAHdb fitting uncertainty
        as the ratio of the residual over the total spectrum area.

        """
        tags = ['err', 'e127', 'e112', 'e77', 'e62', 'e33']

        piecewise = dict.fromkeys(tags)

        range = dict.fromkeys(tags)
        range['err'] = [min(self.grid), max(self.grid)]
        range['e127'] = [754.0, 855.0]
        range['e112'] = [855.0, 1000.0]
        range['e77'] = [1000.0, 1495.0]
        range['e62'] = [1495.0, 1712.0]
        range['e33'] = [2900.0, 3125.0]

        if self.observation:
            for key in piecewise.keys():
                sel = np.where(np.logical_and(
                    self.grid >= range[key][0], self.grid <= range[key][1]))[0]
                total_area = np.trapz(np.array(self.observation)[
                                      sel], x=self.grid[sel])
                if total_area == 0:
                    continue
                resid_area = np.trapz(np.abs(np.array(self.observation)[sel] - self.getfit()[sel]),
                                      x=self.grid[sel])
                piecewise[key] = resid_area / total_area
            return piecewise
