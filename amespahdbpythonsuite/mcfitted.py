#!/usr/bin/env python3

from __future__ import annotations
from typing import Optional

import os
import numpy as np
import matplotlib.pyplot as plt  # type: ignore

from scipy import stats  # type: ignore
from specutils import Spectrum1D  # type: ignore

from amespahdbpythonsuite.amespahdb import AmesPAHdb

message = AmesPAHdb.message


class MCFitted:
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
        self.observation = keywords.get("observation", "")

        if isinstance(d, dict):
            if d.get("type", "") == self.__class__.__name__:
                if "mcfits" not in keywords:
                    self.mcfits = d["mcfits"]
                if "distribution" not in keywords:
                    self.distribution = d["distribution"]
                if "observation" not in keywords:
                    self.observation = d["observation"]

    def get(self) -> dict:
        """
        Calls class: :class:`amespahdbpythonsuite.transitions.Transitions.get`.
        Assigns class variables from inherited dictionary.

        """
        d = {}
        d["type"] = self.__class__.__name__
        d["mcfits"] = self.mcfits
        d["distribution"] = self.distribution
        d["observation"] = self.observation

        return d

    def _getstats(self, d=dict()) -> dict:
        """
        Get statistics for the mcfitted spectra.

        Returns:
            stat : dictionary

        """
        s = stats.describe(d)
        stat = {
            "mean": s.mean,
            "std": np.sqrt(s.variance),
            "skew": s.skewness,
            "kurt": s.kurtosis,
        }

        return stat

    def getfit(self) -> dict:
        """
        Retrieves the mean, std, skewness, and kurtosis spectra.

        """
        mcfits = list()
        for mcfit in self.mcfits:
            mcfits.append(mcfit.getfit())

        return self._getstats(mcfits)

    def getbreakdown(self) -> dict:
        """
        Retrieves the breakdown of the MC fitted PAHs.

        """
        keys = [
            "solo",
            "duo",
            "trio",
            "quartet",
            "quintet",
            "neutral",
            "cation",
            "anion",
            "small",
            "large",
            "pure",
            "nitrogen",
            "nc",
        ]

        results: dict = {key: [] for key in keys}

        # Obtain fit breakdown.
        for fit in self.mcfits:
            bd = fit.getbreakdown()
            for key in keys:
                results[key].append(bd[key])

        ret = dict()
        for key in keys:
            ret[key] = self._getstats(results[key])

        return ret

    def getclasses(self) -> dict:
        """
        Retrieves the spectra of the different classes of the MC fitted PAHs.

        """
        mcclasses = list()
        for mcfit in self.mcfits:
            mcclasses.append(mcfit.getclasses())

        classkeys = list(mcclasses[0].keys())
        _lst_classes: dict = {key: [] for key in classkeys}

        # Create dictionary of lists for each breakdown class.
        for mc in mcclasses:
            for key in classkeys:
                if key in mc.keys():
                    _lst_classes[key].append(mc[key])

        # Create the classes dictionary of dictionaries
        classes = dict()

        for key in classkeys:
            classes[key] = self._getstats(_lst_classes[key])

        return classes

    def plot(self, **keywords):
        """
        Plot the MC sampled fit and breakdown components.

        """
        from astropy.nddata import StdDevUncertainty  # type: ignore

        if keywords.get("datalabel", False):
            datalabel = keywords["datalabel"]
        else:
            datalabel = "obs"

        obs = self.getobservation()

        # Get MC average fit and breakdown spectra.
        fit = self.getfit()
        components = self.getclasses()

        # Plot.
        fig, ax = plt.subplots()

        if keywords.get("wavelength", False):
            x = 1e4 / obs.spectral_axis.value
            xtitle = "Wavelength [micron]"
        else:
            x = obs.spectral_axis.value
            xtitle = (
                self.mcfits[0].units["abscissa"]["label"]
                + " ["
                + self.mcfits[0].units["abscissa"]["unit"].to_string("latex_inline")
                + "]"
            )

        ax.minorticks_on()
        ax.tick_params(which="major", right="on", top="on", direction="in", length=5)
        ax.tick_params(which="minor", right="on", top="on", direction="in", length=3)
        ax.set_xlim((min(x), max(x)))
        ax.set_xlabel(f"{xtitle}")
        ax.set_ylabel(
            self.mcfits[0].units["ordinate"]["label"]
            + " ["
            + self.mcfits[0].units["ordinate"]["unit"].to_string("latex_inline")
            + "]",
        )

        if isinstance(obs.uncertainty, StdDevUncertainty):
            ax.errorbar(
                x,
                obs.flux.value,
                yerr=obs.uncertainty.array,
                fmt="o",
                mfc="white",
                color="k",
                ecolor="k",
                markersize=3,
                elinewidth=0.2,
                capsize=0.8,
                label=datalabel,
            )

        else:
            ax.scatter(x, obs.flux.value, color="k", s=5, label=datalabel)

        ax.plot(x, fit["mean"], color="tab:purple", linewidth=1.2, label="fit")
        ax.fill_between(
            x,
            fit["mean"] - fit["std"],
            fit["mean"] + fit["std"],
            color="tab:purple",
            alpha=0.3,
        )

        if keywords.get("charge"):
            ptype = "charge"
            if isinstance(components["anion"]["mean"], np.ndarray):
                ax.plot(
                    x,
                    components["anion"]["mean"],
                    color="tab:red",
                    linewidth=1.2,
                    label="anion",
                )
                ax.fill_between(
                    x,
                    components["anion"]["mean"] - components["anion"]["std"],
                    components["anion"]["mean"] + components["anion"]["std"],
                    color="tab:red",
                    alpha=0.3,
                )

            if isinstance(components["neutral"]["mean"], np.ndarray):
                ax.plot(
                    x,
                    components["neutral"]["mean"],
                    color="tab:green",
                    linewidth=1.2,
                    label="neutral",
                )
                ax.fill_between(
                    x,
                    components["neutral"]["mean"] - components["neutral"]["std"],
                    components["neutral"]["mean"] + components["neutral"]["std"],
                    color="tab:green",
                    alpha=0.3,
                )

            if isinstance(components["cation"]["mean"], np.ndarray):
                ax.plot(
                    x,
                    components["cation"]["mean"],
                    color="tab:blue",
                    linewidth=1.2,
                    label="cation",
                )
                ax.fill_between(
                    x,
                    components["cation"]["mean"] - components["cation"]["std"],
                    components["cation"]["mean"] + components["cation"]["std"],
                    color="tab:blue",
                    alpha=0.3,
                )

        elif keywords.get("size"):
            ptype = "size"
            if isinstance(components["small"]["mean"], np.ndarray):
                ax.plot(x, components["small"]["mean"], color="tab:red", label="small")
                ax.fill_between(
                    x,
                    components["small"]["mean"] - components["small"]["std"],
                    components["small"]["mean"] + components["small"]["std"],
                    color="tab:red",
                    alpha=0.3,
                )

            if isinstance(components["large"]["mean"], np.ndarray):
                ax.plot(
                    x, components["large"]["mean"], color="tab:green", label="large"
                )
                ax.fill_between(
                    x,
                    components["large"]["mean"] - components["large"]["std"],
                    components["large"]["mean"] + components["large"]["std"],
                    color="tab:green",
                    alpha=0.3,
                )

        elif keywords.get("composition"):
            ptype = "composition"
            if isinstance(components["pure"]["mean"], np.ndarray):
                ax.plot(x, components["pure"]["mean"], color="tab:red", label="pure")
                ax.fill_between(
                    x,
                    components["pure"]["mean"] - components["pure"]["std"],
                    components["pure"]["mean"] + components["pure"]["std"],
                    color="tab:red",
                    alpha=0.3,
                )

            if isinstance(components["nitrogen"]["mean"], np.ndarray):
                ax.plot(
                    x,
                    components["nitrogen"]["mean"],
                    color="tab:green",
                    label="nitrogen",
                )
                ax.fill_between(
                    x,
                    components["nitrogen"]["mean"] - components["nitrogen"]["std"],
                    components["nitrogen"]["mean"] + components["nitrogen"]["std"],
                    color="tab:green",
                    alpha=0.3,
                )

        else:
            ptype = "fitted"

        ax.axhline(0, linestyle="--", color="gray")
        ax.legend(fontsize=10)

        if keywords.get("save", False):
            if keywords.get("output"):
                if os.path.isdir(keywords["output"]):
                    fig.savefig(
                        f"{keywords['output']}/mc_{ptype}_breakdown.{keywords['ftype']}"
                    )
                else:
                    fig.savefig(
                        f"{keywords['output']}_mc_{ptype}_breakdown.{keywords['ftype']}"
                    )
            else:
                fig.savefig(f"mc_{ptype}_breakdown.{keywords['ftype']}")
        else:
            plt.show()

        plt.close(fig)

    def write(self, filename: str = "") -> None:
        """
        Write the spectra to file as an IPAC-table.

        """
        import sys
        import datetime
        from astropy.table import Table  # type: ignore

        if not filename:
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
            "SAMPLES": f"{len(self.mcfits)}",
        }

        for key, value in kv.items():
            if not value.isnumeric():
                hdr.append(f"{key:8} = '{value}'")
            else:
                hdr.append(f"{key:8} = {value}")

        tbl = Table(
            names=("attribute", "mean", "std", "skew", "kurt"),
            dtype=(
                "U25",
                "float64",
                "float64",
                "float64",
                "float64",
            ),
            meta={"comments": hdr},
        )

        for key, vals in self.getbreakdown().items():
            tbl.add_row([key, vals["mean"], vals["std"], vals["skew"], vals["kurt"]])

        tbl.write(filename, format="ipac", overwrite=True)

        message(f"WRITTEN: {filename}")

    def geterror(self) -> dict:
        """
        Obtains the PAHdb fitting uncertainty from the fitted geterror method,
        as the ratio of the residual over the total spectrum area.

        """
        mcerr = list()
        for mcfit in self.mcfits:
            mcerr.append(mcfit.geterror())

        errkeys = list(mcerr[0].keys())
        _lst_classes: dict = {key: [] for key in errkeys}

        # Create dictionary of lists for each error.
        for mc in mcerr:
            for key in errkeys:
                if key in mc.keys():
                    _lst_classes[key].append(mc[key])

        # Create the errors dictionary of dictionaries
        errors: dict = dict()

        for key in errkeys:
            if None in _lst_classes[key]:
                errors[key] = None
            else:
                errors[key] = self._getstats(_lst_classes[key])

        return errors

    def getobservation(self) -> Spectrum1D:
        """
        Retrieves the observation.

        """
        return self.observation

    def getresidual(self) -> dict:
        """
        Retrieves the statistics spectra for the residuals of the MC fits.
        """
        mcres = list()
        for mcfit in self.mcfits:
            mcres.append(mcfit.getresidual())

        return self._getstats(mcres)
