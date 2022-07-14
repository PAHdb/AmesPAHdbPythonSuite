#!/usr/bin/env python3

from __future__ import annotations
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt  # type: ignore

from amespahdbpythonsuite.amespahdb import AmesPAHdb

message = AmesPAHdb.message


class MCfitted:
    """
    AmesPAHdbPythonSuite mcfitted class.
    Contains methods to fit and plot the input spectrum using a Monte Carlo approach.

    """

    def __init__(self, d: Optional[dict] = None, **keywords) -> None:
        self.__set(d, **keywords)

    def set(self, d: Optional[dict] = None, **keywords) -> None:
        """
        Calls class: :class:`amespahdbpythonsuite.transitions.Transitions.set` to parse keywords.

        """
        self.__set(d, **keywords)

    def __set(self, d: Optional[dict] = None, **keywords) -> None:
        """
        Populate data dictionary helper.

        """
        if isinstance(d, dict):
            if d.get("type", "") == self.__class__.__name__:
                if "mcfits" not in keywords:
                    self.mcfits = d["mcfits"]
                if "mcpredicted" not in keywords:
                    self.mcpredicted = d["mcpredicted"]
                if "distribution" not in keywords:
                    self.distribuion = d["distribuion"]

        self.mcfits = keywords.get("mcfits", list())
        self.mcpredicted = keywords.get("mcpredicted", False)
        self.distribution = keywords.get("distribution", "")

    def get(self) -> dict:
        """
        Calls class: :class:`amespahdbpythonsuite.transitions.Transitions.get`.
        Assigns class variables from inherited dictionary.

        """
        d = {}
        d["type"] = self.__class__.__name__
        d["mcfits"] = self.mcfits
        d["mcpredicted"] = self.mcpredicted
        d["distribution"] = self.distribution

        return d

    def getstats(self, save=False, **keywords):
        """
        Obtain statistics (min, max, median, mean, std) for the MC fitted parameters,
        and save to file.

        """

        if 'results' not in keywords:
            results = self.getbreakdown()

        rkeys = list(results.keys())

        if save:
            # Save MC breakdown stats to file.
            o = open('mcfitted_statistics.txt', 'w')
            o.write('# param min max median mean std\n')

            for i, key in enumerate(rkeys):
                o.write(f'{key} {np.min(results[key]):.5f} '
                        f'{np.max(results[key]):.5f} '
                        f'{np.median(results[key]):.5f} '
                        f'{np.mean(results[key]):.5f} '
                        f'{np.std(results[key]):.5f}\n')
            o.close()
        else:
            print('# param min max median mean std')
            for i, key in enumerate(rkeys):
                print(f'{key} {np.min(results[key]):.5f} '
                      f'{np.max(results[key]):.5f} '
                      f'{np.median(results[key]):.5f} '
                      f'{np.mean(results[key]):.5f} '
                      f'{np.std(results[key]):.5f}')

    def averagespec(self, **keywords):
        """
        Calculate average and std spectra.

        """
        #  Getting the average spectra
        components = {}
        lst_classes = {}
        avg_classes = {}
        std_classes = {}

        # Get average and std spectrum.
        if 'fits' not in keywords:
            fits = self.getfit()
        avg_fit = np.mean(fits, axis=0)
        std_fit = np.std(fits, axis=0)

        components['fit'] = {'mean_spec': avg_fit, 'std_spec': std_fit}

        # Get average and std spectrum of classes.
        if 'classes' not in keywords:
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
            std_classes[key] = np.std(lst_classes[key], axis=0)
            components[key] = {'mean_spec': avg_classes[key], 'std_spec': std_classes[key]}

        if 'predicted' in keywords:
            avg_pred = np.mean(self.mcpredicted, axis=0)
            std_pred = np.std(self.mcpredicted, axis=0)
            components['predicted'] = {'mean_spec': avg_pred, 'std_spec': std_pred}

        return components

    def getfit(self, **keywords) -> list:
        """
        Retrieves the sum of fit values.

        """
        fits = []
        for fit in self.mcfits:
            fits.append(fit.getfit())

        return fits

    def getbreakdown(self, **keywords) -> dict:
        """
        Retrieves the breakdown of the MC fitted PAHs.

        """
        results = {'solo': [],
                   'duo': [],
                   'trio': [],
                   'quartet': [],
                   'quintet': [],
                   'anion': [],
                   'neutral': [],
                   'cation': [],
                   'small': [],
                   'large': [],
                   'nitrogen': [],
                   'pure': [],
                   'nc': [],
                   'error': []
                   }

        keys = list(results.keys())

        # Obtain fit breakdown.
        for fit in self.mcfits:
            bd = fit.getbreakdown()
            for key in keys[:-1]:
                results[key].append(bd[key])
            # Obtain fit uncertainty.
            error = fit.geterror()['err']
            results['error'].append(error)

        return results

    def getclasses(self, **keywords) -> list:
        """
        Retrieves the spectra of the different classes of the MC fitted PAHs.

        """
        self.mcclasses = []
        for fit in self.mcfits:
            self.mcclasses.append(fit.getclasses())

        return self.mcclasses

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

        ax.plot(x, components['fit']['mean_spec'], color='tab:purple', linewidth=1.2, label='fit')
        ax.fill_between(x, components['fit']['mean_spec'] - components['fit']['std_spec'],
                        components['fit']['mean_spec'] + components['fit']['std_spec'],
                        color='tab:purple', alpha=0.3)

        if keywords.get('charge'):
            ptype = 'charge'
            if isinstance(components['anion']['mean_spec'], np.ndarray):
                ax.plot(x, components['anion'], color='tab:red', linewidth=1.2, label='anion')
                ax.fill_between(x, components['anion']['mean_spec'] - components['anion']['std_spec'],
                                components['anion']['mean_spec'] + components['anion']['std_spec'],
                                color='tab:red', alpha=0.3)

            if isinstance(components['neutral']['mean_spec'], np.ndarray):
                ax.plot(x, components['neutral']['mean_spec'], color='tab:green', linewidth=1.2, label='neutral')
                ax.fill_between(x, components['neutral']['mean_spec'] - components['neutral']['std_spec'],
                                components['neutral']['mean_spec'] + components['neutral']['std_spec'],
                                color='tab:green', alpha=0.3)

            if isinstance(components['cation']['mean_spec'], np.ndarray):
                ax.plot(x, components['cation'], color='tab:blue', linewidth=1.2, label='cation')
                ax.fill_between(x, components['cation']['mean_spec'] - components['cation']['std_spec'],
                                components['cation']['mean_spec'] + components['cation']['std_spec'],
                                color='tab:blue', alpha=0.3)

        if keywords.get('size'):
            ptype = 'size'
            if isinstance(components['small']['mean_spec'], np.ndarray):
                ax.plot(x, components['small']['mean_spec'], color='tab:red', label='small')
                ax.fill_between(x, components['small']['mean_spec'] - components['small']['std_spec'],
                                components['small']['mean_spec'] + components['small']['std_spec'],
                                color='tab:red', alpha=0.3)

            if isinstance(components['large']['mean_spec'], np.ndarray):
                ax.plot(x, components['large']['mean_spec'], color='tab:green', label='large')
                ax.fill_between(x, components['large']['mean_spec'] - components['large']['std_spec'],
                                components['large']['mean_spec'] + components['large']['std_spec'],
                                color='tab:green', alpha=0.3)

        if keywords.get('composition'):
            ptype = 'composition'
            if isinstance(components['pure']['mean_spec'], np.ndarray):
                ax.plot(x, components['pure']['mean_spec'], color='tab:red', label='pure')
                ax.fill_between(x, components['pure']['mean_spec'] - components['pure']['std_spec'],
                                components['pure']['mean_spec'] + components['pure']['std_spec'],
                                color='tab:red', alpha=0.3)

            if isinstance(components['nitrogen']['mean_spec'], np.ndarray):
                ax.plot(x, components['nitrogen']['mean_spec'], color='tab:green', label='nitrogen')
                ax.fill_between(x, components['nitrogen']['mean_spec'] - components['nitrogen']['std_spec'],
                                components['nitrogen']['mean_spec'] + components['nitrogen']['std_spec'],
                                color='tab:green', alpha=0.3)

        ax.axhline(0, linestyle='--', color='gray')
        ax.legend(fontsize=10)

        if keywords.get('save'):
            if not isinstance(ptype, str):
                ptype = 'fitted'
            fig.savefig(f"mc_{ptype}_breakdown.pdf")
        else:
            plt.show()

        plt.close(fig)
