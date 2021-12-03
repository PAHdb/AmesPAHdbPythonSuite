#!/usr/bin/env python3

import numpy as np
import multiprocessing
import copy
import time
from datetime import timedelta
from functools import partial

from scipy import integrate
from scipy import optimize

from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite.data import Data

message = AmesPAHdb.message


class Transitions(Data):
    """
    AmesPAHdbPythonSuite transitions class.
    Contains methods to create and apply an emission model.

    """

    def __init__(self, d=None, **keywords):
        """
        Initialize transitions class.

        """
        self.__shift = 0.0
        self.set(d, **keywords)

    def set(self, d=None, **keywords):
        """
        Calls class: :class:`amespahdbpythonsuite.data.Data.set to parse keywords.
        Checks type of database object (e.g., theoretical, experimental)

        """
        super().set(d, **keywords)
        if d:
            if d.get('type', '') == self.__class__.__name__:
                if not keywords.get('shift'):
                    self.__shift = d['shift']
                d['type'] = 'Transitions'

        if keywords.get('shift'):
            self.shift = keywords.get('shift')

    def get(self):
        """
        Calls class: :class:`amespahdbpythonsuite.data.Data.get to get keywords.

        """
        d = super().get()
        d['type'] = self.__class__.__name__
        d['shift'] = self.__shift
        return copy.deepcopy(d)

    def shift(self, shift):
        """
        Shifts transitions frequency by provided value.

        """
        self.__shift += shift

        for key in self.data:
            for d in self.data[key]:
                d['frequency'] += shift
        message(f'TOTAL SHIFT: {self.__shift} /cm')

    def fixedtemperature(self, t):
        """
        Applies the Fixed Temperature emission model.

        Parameters
        ----------
        t : float
            Excitation temperature in Kelvin.

        """
        if self.model:
            if self.model['type'] != 'zerokelvin_m':
                message(f'AN EMISSION MODEL HAS ALREADY BEEN APPLIED: {self.model["type"]}')
                return

        message('APPLYING FIXED TEMPERATURE EMISSION MODEL')

        for uid in self.uids:
            f = np.array([d['frequency'] for d in self.data[uid]])

            intensity = 2.4853427121856266e-23 * f ** 3 / (np.exp(1.4387751297850830401 * f / t) - 1.0)
            for (d, i) in zip(self.data[uid], intensity):
                d['intensity'] *= i

    def calculatedtemperature(self, e, **keywords):
        """
        Applies the Calculated Temperature emission model.

        Parameters
        ----------
        e : float
            Excitation energy in erg.

        """
        if self.type != 'theoretical':
            message('THEORETICAL DATABASE REQUIRED FOR EMISSION MODEL')
            return

        if self.model:
            if self.model['type'] != 'zerokelvin_m':
                message(f'AN EMISSION MODEL HAS ALREADY BEEN APPLIED: {self.model["type"]}')
                return

        message('APPLYING CALCULATED TEMPERATURE EMISSION MODEL')

        global energy

        energy = e

        self.model = {'type': 'calculatedtemperature_m',
                      'energy': e,
                      'temperatures': [],
                      'description': ''}

        self.units['ordinate'] = {'unit': 2,
                                  'str': 'integrated spectral radiance [erg/s/PAH]'}

        print(57 * '=')

        i = 0

        nuids = len(self.uids)

        for uid in self.uids:
            # Start timer.
            tstart = time.perf_counter()

            print('SPECIES                          : %d/%d' % (i + 1, nuids))
            print('UID                              : %d' % uid)
            print('MEAN ABSORBED ENERGY             : %f +/- %f eV' % (e / 1.6021765e-12, 0.0))

            global frequencies
            frequencies = np.array([d['frequency'] for d in self.data[uid]])

            Tmax = optimize.brentq(self.attainedtemperature, 2.73, 5000.0)

            print('MAXIMUM ATTAINED TEMPERATURE     : %f Kelvin' % Tmax)

            self.model['temperatures'].append({'uid': uid, 'temperature': Tmax})

            for d in self.data[uid]:
                if d['intensity'] > 0:
                    d['intensity'] *= 2.4853427121856266e-23 * \
                        d['frequency'] ** 3 / (np.exp(1.4387751297850830401 * d['frequency'] / Tmax) - 1.0)

            # Stop timer and calculate elapsed time.
            elapsed = timedelta(seconds=(time.perf_counter() - tstart))
            print(f'Elapsed time: {elapsed}')

            i += 1

        print(57 * '=')

    def _cascade_em_model(self, e, uid):
        """
        A partial method of :meth:`amespahdbpythonsuite.transitions.cascade`
        used when multiprocessing is required.

        Parameters
        ----------

        uid : int
            single UID value.

        Returns
        -------
        temp : dict
            Dictionary of Tmax.
        ud : dict
            Dictionary of calculated intensities for the given UID.

        """
        global energy
        energy = e

        global frequencies
        frequencies = np.array([d['frequency'] for d in self.data[uid]])

        global intensities
        intensities = np.array([d['intensity'] for d in self.data[uid]])

        Tmax = optimize.brentq(self.attainedtemperature, 2.73, 5000.0)

        for d in self.data[uid]:
            if d['intensity'] > 0:
                global frequency
                frequency = d['frequency']
                d['intensity'] *= d['frequency'] ** 3 * integrate.quad(self.featurestrength, 2.73, Tmax)[0]

        temp = {'uid': uid, 'temperature': Tmax}
        ud = {uid: self.data[uid]}
        return ud, temp

    def cascade(self, e, **keywords):
        """
        Applies the Cascade emission model.

        Parameters
        ----------
        e : float
            Excitation energy in erg.

        """
        if self.type != 'theoretical':
            message('THEORETICAL DATABASE REQUIRED FOR EMISSION MODEL')
            return

        if self.model:
            if self.model['type'] != 'zerokelvin_m':
                message(f'AN EMISSION MODEL HAS ALREADY BEEN APPLIED: {self.model["type"]}')
                return

        message('APPLYING CASCADE EMISSION MODEL')

        tstart = time.perf_counter()

        global energy

        energy = e

        self.model = {'type': 'cascade_m',
                      'energy': e,
                      'temperatures': [],
                      'description': ''}

        self.units['ordinate'] = {'unit': 3,
                                  'str': 'integrated radiant energy [erg]'}

        print(57 * '=')

        if keywords.get('multiprocessing'):
            cascade_em_model = partial(self._cascade_em_model, e)
            if keywords.get('ncores'):
                ncores = keywords.get('ncores')
                print(f'Using multiprocessing module with {ncores} cores')
                pool = multiprocessing.Pool(processes=ncores)
                data, temp = zip(*pool.map(cascade_em_model, self.uids))
            else:
                print('Using multiprocessing module with 2 cores')
                pool = multiprocessing.Pool(processes=2)
                data, temp = zip(*pool.map(cascade_em_model, self.uids))
            pool.close()
            pool.join()

            # Re-assign self.data.
            self.data = {}
            for d in list(data):
                for k, v in d.items():
                    self.data[k] = v

            self.model['temperatures'] = temp

        else:
            i = 0
            nuids = len(self.uids)
            for uid in self.uids:

                print('SPECIES                          : %d/%d' % (i + 1, nuids))
                print('UID                              : %d' % uid)
                print('MEAN ABSORBED ENERGY             : %f +/- %f eV' % (e / 1.6021765e-12, 0.0))

                global frequencies
                frequencies = np.array([d['frequency'] for d in self.data[uid]])

                global intensities
                intensities = np.array([d['intensity'] for d in self.data[uid]])

                Tmax = optimize.brentq(self.attainedtemperature, 2.73, 5000.0)

                print('MAXIMUM ATTAINED TEMPERATURE     : %f Kelvin' % Tmax)

                self.model['temperatures'].append({'uid': uid, 'temperature': Tmax})

                for d in self.data[uid]:
                    if d['intensity'] > 0:
                        global frequency
                        frequency = d['frequency']
                        d['intensity'] *= d['frequency'] ** 3 * integrate.quad(self.featurestrength, 2.73, Tmax)[0]

                i += 1

            print(57 * '=')

        elapsed = timedelta(seconds=(time.perf_counter() - tstart))
        print(f'Elapsed time: {elapsed}\n')

    def _get_intensities(self, npoints, xmin, xmax, clip, width, x, gaussian, drude, uid):
        """
        A partial method of :meth:`amespahdbpythonsuite.transitions.convolve`
        used when multiprocessing is required.

        Parameters
        ----------

        npoints : int
            Number of grid points.
        xmin : float
            Minimum value of grid.
        xmax : float
            Maximum value of grid.
        clip : float
            Value to clip and define the frequency range of the profile calculation.
        width : float
            Width of the line profile.
        x : array
            Grid array.
        gaussian : str
            String to indicate Gaussian profile
        drude : str
            String to indicate Drude profile
        uid : int
            Single UID value

        Returns
        -------
        ud : dict
            Dictionary of convolutged intensities for the given UID.

        """
        s = np.zeros(npoints)
        f = [v for v in self.data[uid] if v['frequency'] >= xmin - clip * width
             and v['frequency'] <= xmax + clip * width]
        for t in f:
            if t['intensity'] > 0:
                s += t['intensity'] * \
                    self.__lineprofile(x,
                                       t['frequency'],
                                       width,
                                       gaussian=gaussian,
                                       drude=drude)

        ud = {uid: s}
        return ud

    def convolve(self, **keywords):
        """
        Convolve transitions with a line profile.
        Calls class:
        :class:`amespahdbpythonsuite.spectrum.Spectrum` to retrieve the respective object.

        """

        fwhm = keywords.get('fwhm', 15.0)

        if keywords.get('gaussian'):
            width = 0.5 * fwhm / np.sqrt(2.0 * np.log(2.0))
            clip = 3.0
            profile = 'Gaussian'
            message('USING GAUSSIAN LINE PROFILES')

        elif keywords.get('drude'):
            width = 1.0 / fwhm
            clip = 11.0
            profile = 'Drude'
            message('USING DRUDE LINE PROFILES')

        else:
            width = 0.5 * fwhm
            clip = 22.0
            profile = 'Lorentzian'
            message('USING LORENTZIAN LINE PROFILES')

        x = keywords.get('grid', [])
        if not len(x):
            r = keywords.get('xrange', [])
            if len(r):
                xmin = min(r)
                xmax = max(r)
            else:
                xmin = 1.0
                xmax = 4000.0

            npoints = keywords.get('npoints', False)
            if not npoints:
                npoints = 400

            x = np.arange(xmin, xmax, (xmax - xmin) / npoints)

        else:
            x = np.asarray(x)
            xmin = min(x)
            xmax = max(x)
            npoints = len(x)

        message(f'GRID: (XMIN,XMAX)=({xmin:.3f}, {xmax:.3f}); {npoints} POINTS')
        message(f'FWHM: {fwhm} /cm')

        print('Convolving')
        # Start timer.
        tstart = time.perf_counter()

        gaussian = keywords.get('gaussian', False)
        drude = keywords.get('drude', False)
        d = {}

        if keywords.get('multiprocessing'):
            get_intensities = partial(self._get_intensities, npoints,
                                      xmin, xmax, clip, width,
                                      x, gaussian, drude)
            if keywords.get('ncores'):
                ncores = keywords.get('ncores')
                print(f'Using multiprocessing module with {ncores} cores')
                pool = multiprocessing.Pool(processes=ncores)
                data = pool.map(get_intensities, self.uids)
            else:
                print('Using multiprocessing module with 2 cores')
                pool = multiprocessing.Pool(processes=2)
                data = pool.map(get_intensities, self.uids)
            pool.close()
            pool.join()

            # Re-assign self.data.
            for ud in list(data):
                for k, v in ud.items():
                    d[k] = v

        else:
            for uid in self.uids:
                s = np.zeros(npoints)
                f = [v for v in self.data[uid] if
                     (v['frequency'] >= xmin - clip * width) and
                     (v['frequency'] <= xmax + clip * width)]
                for t in f:
                    if t['intensity'] > 0:
                        s += t['intensity'] * self.__lineprofile(x,
                                                                 t['frequency'],
                                                                 width,
                                                                 gaussian=keywords.get('gaussian', False),
                                                                 drude=keywords.get('drude', False))
                d[uid] = s

        elapsed = timedelta(seconds=(time.perf_counter() - tstart))
        print(f'Elapsed time: {elapsed}')

        from amespahdbpythonsuite.spectrum import Spectrum

        if self.model['type'] == 'zerokelvin_m':
            self.units['ordinate'] = {'unit': 2,
                                      'str': 'cross-section [x10$^{5}$ cm$^{2}$/mol]'}
        elif self.model['type'] == 'cascade_m':
            self.units['ordinate'] = {'unit': 3,
                                      'str': 'radiant energy [x10$^{5}$ erg/cm$^{-1}$/mol]'}

        return Spectrum(type=self.type,
                        version=self.version,
                        data=d,
                        pahdb=self.pahdb,
                        uids=self.uids,
                        model=self.model,
                        units={'abscissa': self.units['abscissa'],
                               'ordinate': self.units['ordinate']},
                        shift=self.__shift,
                        grid=x,
                        profile=profile,
                        fwhm=fwhm)

    def __lineprofile(self, x, x0, width, **keywords):
        """
        Calculate Gaussian, Drude, or Lorentzian line profiles.

        Parameters
        ----------
        x : array
            Grid array.
        x0 : float
            Central frequency
        width : float
            Width of the line profile.

        """
        if keywords.get('gaussian', False):
            return (1.0 / (width * np.sqrt(2.0 * np.pi))) * \
                np.exp(-(x - x0) ** 2 / (2.0 * width ** 2))
        elif keywords.get('drude', False):
            return (2.0 / (np.pi * x0 * width)) * \
                width ** 2 / ((x / x0 - x0 / x) ** 2 + width ** 2)
        elif keywords.get('lorentzian', True):
            return (width / np.pi) / ((x - x0) ** 2 + width ** 2)

    def plot(self, **keywords):
        """
        Plot the transitions absorption spectrum.

        """
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm

        _, ax = plt.subplots()
        ax.minorticks_on()
        ax.tick_params(which='major', right='on', top='on', direction='in', length=5)
        ax.tick_params(which='minor', right='on', top='on', direction='in', length=3)
        colors = cm.rainbow(np.linspace(0, 1, len(self.uids)))
        for uid, col in zip(self.uids, colors):
            f = [v for v in self.data[uid]]
            x = [d['frequency'] for d in f]
            y = [d['intensity'] for d in f]
            ax.bar(x, y, 20, color=col, alpha=0.5)
            ax.tick_params(axis='both', labelsize=14)

        plt.gca().invert_xaxis()
        plt.xlabel(self.units['abscissa']['str'], fontsize=14)
        plt.ylabel(self.units['ordinate']['str'], fontsize=14)

        if keywords.get('show'):
            plt.show()
        elif keywords.get('outfile'):
            outfile = keywords.get('outfile')
            plt.savefig(f'{outfile}.pdf', bbox_inches='tight')

    @staticmethod
    def featurestrength(T):
        """
        Calculate a feature's strength covolved with a blackbody.

        Parameters
        ----------
        T : float
            Excitation temperature in Kelvin.

        """
        global frequency
        global frequencies
        global intensities

        val1 = 1.4387751297850830401 * frequency / T

        val2 = 1.4387751297850830401 * frequencies / T

        # TODO: replace np.exp with np.expm1, and avoid overflow.
        return (Transitions.heatcapacity(T) / (np.exp(val1) - 1)) * \
            (1 / np.sum(intensities * (frequencies) ** 3 / (np.exp(val2) - 1)))

    @staticmethod
    def attainedtemperature(T):
        """
        Calculate a PAH's temperature after absorbing a given amount of energy.

        Parameters
        ----------
        T : float
            Excitation temperature in Kelvin.

        """
        global energy

        return integrate.quad(Transitions.heatcapacity, 2.73, T)[0] - energy

    def heatcapacity(T):
        """
        Calculate heat capacity.

        Parameters
        ----------
        T : float
            Excitation temperature in Kelvin.

        """

        global frequencies

        val = 1.4387751297850830401 * frequencies / T

        return 1.3806505e-16 * np.sum(np.exp(-val) * (val / (1.0 - np.exp(-val))) ** 2)
