#!/usr/bin/env python3

from __future__ import annotations
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt  # type: ignore

from scipy import stats

from amespahdbpythonsuite.amespahdb import AmesPAHdb

message = AmesPAHdb.message


class MCfitted:
    """
    AmesPAHdbPythonSuite mcfitted class.
    Contains methods to fit and plot the input spectrum using a Monte Carlo approach.

    """

    def __init__(self, d: Optional[dict] = None, **keywords) -> None:
        self.__set(d, **keywords)

    def __set(self, d: Optional[dict] = None, **keywords) -> None:
        """
        Populate data dictionary helper.

        """
        self.mcfits = keywords.get("mcfits", list())
        self.distribution = keywords.get("distribution", "")

        if isinstance(d, dict):
            if d.get("type", "") == self.__class__.__name__:
                if "mcfits" not in keywords:
                    self.mcfits = d["mcfits"]
                if "distribution" not in keywords:
                    self.distribuion = d["distribuion"]

    def get(self) -> dict:
        """
        Calls class: :class:`amespahdbpythonsuite.transitions.Transitions.get`.
        Assigns class variables from inherited dictionary.

        """
        d = {}
        d["type"] = self.__class__.__name__
        d["mcfits"] = self.mcfits
        d["distribution"] = self.distribution

        return d

    def averagespec(self, **keywords):
        """
        Calculate average and std spectra.

        """
        #  Getting the average spectra
        components = dict()
        lst_classes = dict()
        avg_classes = dict()
        std_classes = dict()

        # Get average and std spectrum.
        components['fit'] = self.getfit()

        # Get average and std spectrum of classes.
        classes = self.getclasses()
        classkeys = list(classes[0].keys())

        for key in classkeys:
            lst_classes[key] = []
            avg_classes[key] = []
            std_classes[key] = []
        for c in classes:
            for key in classkeys:
                if key in c.keys():
                    lst_classes[key].append(c[key])
        for key in classkeys:
            avg_classes[key] = np.mean(lst_classes[key], axis=0)
            std_classes[key] = np.var(lst_classes[key], axis=0)
            components[key] = {'mean': avg_classes[key], 'var': std_classes[key]}

        return components

    def getfit(self, **keywords) -> dict:
        """
        Retrieves the sum of fit values.

        """
        mcfits = list()
        for mcfit in self.mcfits:
            mcfits.append(mcfit.getfit())

        return {'mean' : np.mean(mcfits, axis=0), 'var': np.var(mcfits, axis=0)}

    def getbreakdown(self, **keywords) -> dict:
        """
        Retrieves the breakdown of the MC fitted PAHs.

        """
        keys = ['solo', 'duo', 'trio', 'quartet', 'quintet',
                'neutral', 'cation', 'anion',
                'small', 'large',
                'pure', 'nitrogen',
                'nc', 'error']

        results = {key: [] for key in keys}

        # Obtain fit breakdown.
        for fit in self.mcfits:
            bd = fit.getbreakdown()
            for key in keys[:-1]:
                results[key].append(bd[key])
            # Obtain fit uncertainty.
            error = fit.geterror()['err']
            results['error'].append(error)

        ret = dict()
        stat = dict()
        for k in keys:
            ret[k] = {'mean': np.mean(results[k]), 'var': np.var(results[k])}
            stat[k] = stats.describe(results[k])

        for key, value in stat.items():
            print(key, value)
        
        if keywords.get('write'):
            self.write(stat, keywords.get('write'))

        return ret

    def getclasses(self, **keywords) -> list:
        """
        Retrieves the spectra of the different classes of the MC fitted PAHs.

        """
        mcclasses = list()
        for mcfit in self.mcfits:
            mcclasses.append(mcfit.getclasses())

        return mcclasses

    def plot(self, **keywords):
        """
        Plot the MC sampled fit and breakdown components.

        """
        from astropy.nddata import StdDevUncertainty  # type: ignore

        obs = self.mcfits[0].observation
        units = self.mcfits[0].units

        if 'components' not in keywords:
            components = self.averagespec()

        # Charge
        fig, ax = plt.subplots()

        if keywords.get('wavelength', False):
            x = 1e4 / obs.spectral_axis.value
            xtitle = "Wavelength [micron]"
        else:
            x = obs.spectral_axis.value
            xtitle = (
                units["abscissa"]["label"]
                + " ["
                + units["abscissa"]["unit"].to_string("latex_inline")
                + "]"
            )

        ax.minorticks_on()
        ax.tick_params(which='major', right='on', top='on', direction='in', length=5)
        ax.tick_params(which='minor', right='on', top='on', direction='in', length=3)
        ax.set_xlim((min(x), max(x)))
        ax.set_xlabel(f"{xtitle}")
        ax.set_ylabel(
            units["ordinate"]["label"]
            + " ["
            + units["ordinate"]["unit"].to_string("latex_inline")
            + "]",
        )

        if isinstance(obs.uncertainty, StdDevUncertainty):
            ax.errorbar(x,
                        obs.flux.value,
                        yerr=obs.uncertainty.array,
                        fmt="o",
                        mfc="white",
                        color="k",
                        ecolor="k",
                        markersize=3,
                        elinewidth=0.2,
                        capsize=0.8,
                        label="obs"
                        )

        else:
            ax.scatter(x,
                       obs.flux.value,
                       color='k',
                       s=10,
                       label='obs'
                       )

        ax.plot(x, components['fit']['mean'], color='tab:purple', linewidth=1.2, label='fit')
        ax.fill_between(x, components['fit']['mean'] - components['fit']['var'],
                        components['fit']['mean'] + components['fit']['var'],
                        color='tab:purple', alpha=0.3)

        if keywords.get('charge'):
            ptype = 'charge'
            if isinstance(components['anion']['mean'], np.ndarray):
                ax.plot(x, components['anion'], color='tab:red', linewidth=1.2, label='anion')
                ax.fill_between(x, components['anion']['mean'] - components['anion']['var'],
                                components['anion']['mean'] + components['anion']['var'],
                                color='tab:red', alpha=0.3)

            if isinstance(components['neutral']['mean'], np.ndarray):
                ax.plot(x, components['neutral']['mean'], color='tab:green', linewidth=1.2, label='neutral')
                ax.fill_between(x, components['neutral']['mean'] - components['neutral']['var'],
                                components['neutral']['mean'] + components['neutral']['var'],
                                color='tab:green', alpha=0.3)

            if isinstance(components['cation']['mean'], np.ndarray):
                ax.plot(x, components['cation'], color='tab:blue', linewidth=1.2, label='cation')
                ax.fill_between(x, components['cation']['mean'] - components['cation']['var'],
                                components['cation']['mean'] + components['cation']['var'],
                                color='tab:blue', alpha=0.3)

        if keywords.get('size'):
            ptype = 'size'
            if isinstance(components['small']['mean'], np.ndarray):
                ax.plot(x, components['small']['mean'], color='tab:red', label='small')
                ax.fill_between(x, components['small']['mean'] - components['small']['var'],
                                components['small']['mean'] + components['small']['var'],
                                color='tab:red', alpha=0.3)

            if isinstance(components['large']['mean'], np.ndarray):
                ax.plot(x, components['large']['mean'], color='tab:green', label='large')
                ax.fill_between(x, components['large']['mean'] - components['large']['var'],
                                components['large']['mean'] + components['large']['var'],
                                color='tab:green', alpha=0.3)

        if keywords.get('composition'):
            ptype = 'composition'
            if isinstance(components['pure']['mean'], np.ndarray):
                ax.plot(x, components['pure']['mean'], color='tab:red', label='pure')
                ax.fill_between(x, components['pure']['mean'] - components['pure']['var'],
                                components['pure']['mean'] + components['pure']['var'],
                                color='tab:red', alpha=0.3)

            if isinstance(components['nitrogen']['mean'], np.ndarray):
                ax.plot(x, components['nitrogen']['mean'], color='tab:green', label='nitrogen')
                ax.fill_between(x, components['nitrogen']['mean'] - components['nitrogen']['var'],
                                components['nitrogen']['mean'] + components['nitrogen']['var'],
                                color='tab:green', alpha=0.3)

        ax.axhline(0, linestyle='--', color='gray')
        ax.legend(fontsize=10)

        if keywords.get('save'):
            if not isinstance(ptype, str):
                ptype = 'fitted'
            fig.savefig(f"{keywords['save']}mc_{ptype}_breakdown.pdf")
        else:
            plt.show()

        plt.close(fig)

    def write(self, stat=dict, filename: str = "") -> None:
        """
        Write the spectra to file as an IPAC-table.

        """
        import sys
        import datetime
        from astropy.table import Table  # type: ignore

        if filename == "":
            filename = self.__class__.__name__ + ".tbl"

        hdr = list()

        kv = {
            "DATE": datetime.datetime.now()
            .astimezone()
            .replace(microsecond=0)
            .isoformat(),
            "ORIGIN": "NASA Ames Research Center",
            "CREATOR": f"Python {sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
            "SOFTWARE": "AmesPAHdbPythonSuite",
            "AUTHOR": "Dr. C. Boersma",
            "TYPE": self.__class__.__name__.upper(),
        }

        for key, value in kv.items():
            if not value.isnumeric():
                hdr.append(f"{key:8} = '{value}'")
            else:
                hdr.append(f"{key:8} = {value}")

        tbl = Table(names=('attribute', 'min', 'max', 'mean', 'var', 'skew', 'kurt'),
                    dtype=('U25', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64',),
                    meta={"comments": hdr})

        for key, vals in stat.items():
            tbl.add_row([key, vals[1][0], vals[1][1], vals[2], vals[3], vals[4], vals[5]])

        tbl.write(filename, format='ipac', overwrite=True)

        message(f"WRITTEN: {filename}")
