#!/usr/bin/env python3

from __future__ import annotations
from typing import TYPE_CHECKING, Optional, Union

import numpy as np

from specutils import Spectrum1D, manipulation  # type: ignore
from astropy.nddata import StdDevUncertainty  # type: ignore
import astropy.units as u  # type: ignore
from scipy import optimize  # type: ignore

from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite.transitions import Transitions

if TYPE_CHECKING:
    from amespahdbpythonsuite.coadded import Coadded
    from amespahdbpythonsuite.fitted import Fitted


message = AmesPAHdb.message


class Mcfitted(Transitions):
    """
    AmesPAHdbPythonSuite mcfitted class.
    Contains methods to fit and plot the input spectrum using a Monte Carlo approach.

    """

    def __init__(self, d: Optional[dict] = None, **keywords) -> None:
        super().__init__(d, **keywords)
        self.__set(d, **keywords)

    def set(self, d: Optional[dict] = None, **keywords) -> None:
        """
        Calls class: :class:`amespahdbpythonsuite.transitions.Transitions.set` to parse keywords.

        """
        Transitions.set(self, d, **keywords)
        self.__set(d, **keywords)

    def __set(self, d: Optional[dict] = None, **keywords) -> None:
        """
        Populate data dictionary helper.

        """
        if isinstance(d, dict):
            if d.get("type", "") == self.__class__.__name__:
                if "grid" not in keywords:
                    self.grid = d["grid"]
                if "profile" not in keywords:
                    self.profile = d["profile"]
                if "fwhm" not in keywords:
                    self.fwhm = d["fwhm"]

        self.grid = keywords.get("grid", list())
        self.profile = keywords.get("profile", "")
        self.fwhm = keywords.get("fwhm", 0.0)

    def get(self) -> dict:
        """
        Calls class: :class:`amespahdbpythonsuite.transitions.Transitions.get`.
        Assigns class variables from inherited dictionary.

        """
        d = Transitions.get(self)
        d["type"] = self.__class__.__name__
        d["grid"] = self.grid
        d["profile"] = self.profile
        d["fwhm"] = self.fwhm

        return d

    def stats(self, results):
        """
        Obtain statistics (min, max, median, mean, std) for the MC fitted parameters.

        """

        # Save MC breakdown stats to file.
        o = open('mcfitted_statistics.txt', 'w')
        o.write('# param min max median mean std\n')

        keys = list(results.keys())
        for i, key in enumerate(keys):
            if i < 5:
                o.write(f'{key} {np.min(results[key]):d} {np.max(results[key]):d} {np.median(results[key]):.5f} {np.mean(results[key]):.5f} {np.std(results[key]):.5f}\n')
            else:
                o.write(f'{key} {np.min(results[key]):.5f} {np.max(results[key]):.5f} {np.median(results[key]):.5f} {np.mean(results[key]):.5f} {np.std(results[key]):.5f}\n')
        o.close()

    def averagespec(self, fits, classes, weights):
        import pickle
        #  Getting the average spectra
        components = {}
        lst_classes = {}
        avg_classes = {}
        std_classes = {}

        # Average and std spectrum of fitted and predicted spectra.
        avg_fit = np.mean(fits, axis=0)
        std_fit = np.std(fits, axis=0)

        components['fit'] = {'mean_spec': avg_fit, 'std_spec': std_fit}

        classkeys = ['anion', 'neutral', 'cation', 'small', 'large', 'nitrogen', 'pure']

        for key in classkeys:
            lst_classes[key] = []
            avg_classes[key] = []
            std_classes[key] = []
        for c in classes:
            for key in classkeys:
                lst_classes[key].append(c[key])
        for key in classkeys:
            avg_classes[key] = np.mean(lst_classes[key], axis=0)
            std_classes[key] = np.std(lst_classes[key], axis=0)
            components[key] = {'mean_spec': avg_classes[key], 'std_spec': std_classes[key]}

        # Dump average components pickle
        pickle.dump(components, open('average_mc_spectra.pkl', 'wb'))

        # Calculate average weights.
        luids = [d['uid'] for d in weights]
        unq_uids = np.unique([item for sublist in luids for item in sublist])
        mc_uids = []
        avg_fweight = []
        std_fweight = []

        for uid in unq_uids:
            mc_uids.append(uid)
            w_tmp = []
            for sub in weights:
                idx = np.where(sub['uid'] == uid)[0]
                if len(idx) != 0:
                    w_tmp.append(sub['weights'][idx][0])
                else:
                    w_tmp.append(0.0)
            avg_fweight.append(np.average(w_tmp))
            std_fweight.append(np.std(w_tmp))

        ascii.write([mc_uids, avg_fweight, std_fweight],
                    f'{self.pahdbdir}/{self.galaxy}_{spectype}_mc_avg_fweights.txt',
                    names=['#UID', 'avg_fweight', 'std_fweight'],
                    overwrite=True)
