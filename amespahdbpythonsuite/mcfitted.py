#!/usr/bin/env python3

from __future__ import annotations
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt

from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite.spectrum import Spectrum

message = AmesPAHdb.message


class Mcfitted(Spectrum):
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
        Spectrum.set(self, d, **keywords)
        self.__set(d, **keywords)

    def __set(self, d: Optional[dict] = None, **keywords) -> None:
        """
        Populate data dictionary helper.

        """
        if isinstance(d, dict):
            if d.get("type", "") == self.__class__.__name__:
                if "observation" not in keywords:
                    self.observation = d["observation"]
                if "mcfits" not in keywords:
                    self.mcfits = d["mcfits"]
                if "mcweights" not in keywords:
                    self.mcweights = d["mcweights"]
                if "mcclasses" not in keywords:
                    self.mcclasses = d["mcclasses"]
                if "mcresults" not in keywords:
                    self.mcresults = d["mcresults"]
                if "mcpredicted" not in keywords:
                    self.mcpredicted = d["mcpredicted"]

        self.observation = keywords.get("observation", None)
        self.mcfits = keywords.get("mcfits", list())
        self.mcweights = keywords.get("mcweights", list())
        self.mcclasses = keywords.get("mcclasses", list())
        self.mcresults = keywords.get("mcresults", list())
        self.mcpredicted = keywords.get("mcpredicted", False)

    def get(self) -> dict:
        """
        Calls class: :class:`amespahdbpythonsuite.transitions.Transitions.get`.
        Assigns class variables from inherited dictionary.

        """
        d = Spectrum.get(self)
        d["type"] = self.__class__.__name__
        d["observation"] = self.observation
        d["mcfits"] = self.mcfits
        d["mcweights"] = self.mcweights
        d["mcclasses"] = self.mcclasses
        d["mcresults"] = self.mcresults
        d["mcpredicted"] = self.mcpredicted

        return d

    def stats(self, results):
        """
        Obtain statistics (min, max, median, mean, std) for the MC fitted parameters,
        and save to file.

        """

        # Save MC breakdown stats to file.
        o = open('mcfitted_statistics.txt', 'w')
        o.write('# param min max median mean std\n')

        rkeys = list(results.keys())
        for i, key in enumerate(rkeys):
            o.write(f'{key} {np.min(results[key]):.5f} '
                    f'{np.max(results[key]):.5f} '
                    f'{np.median(results[key]):.5f} '
                    f'{np.mean(results[key]):.5f} '
                    f'{np.std(results[key]):.5f}\n')
        o.close()

    def averagespec(self, mcfits, mcclasses, **keywords):
        """
        Calculate average and std spectra.

        """
        #  Getting the average spectra
        components = {}
        lst_classes = {}
        avg_classes = {}
        std_classes = {}

        # Average and std spectrum.
        avg_fit = np.mean(mcfits, axis=0)
        std_fit = np.std(mcfits, axis=0)

        if 'mcpredicted' in keywords:
            avg_pred = np.mean(self.mcpredicted, axis=0)
            std_pred = np.std(self.mcpredicted, axis=0)
            components['predicted'] = {'mean_spec': avg_pred, 'std_spec': std_pred}

        components['fit'] = {'mean_spec': avg_fit, 'std_spec': std_fit}

        classkeys = ['anion', 'neutral', 'cation', 'small', 'large', 'nitrogen', 'pure']

        for key in classkeys:
            lst_classes[key] = []
            avg_classes[key] = []
            std_classes[key] = []
        for c in mcclasses:
            for key in classkeys:
                if key in c.keys():
                    lst_classes[key].append(c[key])
        for key in classkeys:
            avg_classes[key] = np.mean(lst_classes[key], axis=0)
            std_classes[key] = np.std(lst_classes[key], axis=0)
            components[key] = {'mean_spec': avg_classes[key], 'std_spec': std_classes[key]}

        # Dump average components pickle
        if 'dump' in keywords:
            import pickle
            pickle.dump(components, open('average_mc_spectra.pkl', 'wb'))

        return components

    def plot(self, components, **keywords):
        """
        Plot the MC sampled fit and breakdown components.

        """

        # Charge
        fig, ax = plt.subplots()

        if 'wavelength' in keywords:
            x = 1e4 / self.grid
            ax.set_xlim((min(x), max(x)))
            xtitle = "Wavelength"
        else:
            x = self.grid
            ax.set_xlim((max(x), min(x)))
            xtitle = (
                self.units["abscissa"]["label"]
                + " ["
                + self.units["abscissa"]["unit"].to_string("latex_inline")
                + "]"
            )

        ax.minorticks_on()
        ax.tick_params(which='major', right='on', top='on', direction='in', length=5)
        ax.tick_params(which='minor', right='on', top='on', direction='in', length=3)
        ax.set_xlim((min(x), max(x)))
        ax.set_xlabel(f"{xtitle}")
        ax.set_ylabel(
            self.units["ordinate"]["label"]
            + " ["
            + self.units["ordinate"]["unit"].to_string("latex_inline")
            + "]",
        )

        ax.errorbar(
            x,
            self.observation.spectrum.flux.value,
            yerr=keywords["sigma"],
            fmt="o",
            mfc="white",
            color="k",
            ecolor="k",
            markersize=3,
            elinewidth=0.2,
            capsize=0.8,
            label="obs",
        )

        ax.plot(x, components['fit']['mean_spec'], color='tab:purple', linewidth=1.2, label='fit')
        ax.fill_between(x, components['fit']['mean_spec'] - components['fit']['std_spec'],
                        components['fit']['mean_spec'] + components['fit']['std_spec'],
                        color='tab:purple', alpha=0.3)

        if 'charge' in keywords:
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

        if 'size' in keywords:
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

        if 'composition' in keywords:
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
        ax.legend(fontsize=12)

        basename = keywords['ptype']
        if basename:
            if not isinstance(basename, str):
                basename = 'fitted'
            fig.savefig(f"mc_{keywords['ptype']}_breakdown.{keywords['ftype']}")
        elif 'show' in keywords:
            plt.show()

        plt.close(fig)
