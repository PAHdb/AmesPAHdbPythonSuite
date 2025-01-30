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

        self._fit: Optional[dict] = None
        self._breakdown: Optional[dict] = None
        self._classes: Optional[dict] = None
        self._error: Optional[dict] = None
        self._residual: Optional[dict] = None

    def get(self) -> dict:
        """
        Calls class: :class:`amespahdbpythonsuite.transitions.Transitions.get`.
        Returns a dictionary including PAH UIDs and their corresponding weights.
        """

        return {
            "type": self.__class__.__name__,
            "mcfits": [
                {
                    "uids": fitted.uids,  # List of PAH UIDs in the fit
                    "weights": fitted.weights,  # Dictionary {uid: weight}
                }
                for fitted in self.mcfits
            ],
            "distribution": self.distribution,
            "observation": self.observation,
        }


    def _getstats(self, d=list()) -> dict:
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
        if self._fit is None:
            self._fit = self._getstats([mcfit.getfit() for mcfit in self.mcfits])

        return self._fit

    def getbreakdown(self) -> dict:
        """
        Retrieves the breakdown of the MC fitted PAHs.

        """
        if self._breakdown is None:
            mcfits = iter(self.mcfits)
            mcfit = next(mcfits)
            breakdown = mcfit.getbreakdown()
            results: dict = {k: [] for k in breakdown.keys()}
            for key, val in breakdown.items():
                results[key].append(val)
            for mcfit in mcfits:
                breakdown = mcfit.getbreakdown()
                for key, val in breakdown.items():
                    results[key].append(val)
            self._breakdown = {key: self._getstats(val) for key, val in results.items()}

        return self._breakdown

    def getclasses(self) -> dict:
        """
        Retrieves the spectra of the different classes of the MC fitted PAHs.

        """
        if self._classes is None:
            mcfits = iter(self.mcfits)
            mcfit = next(mcfits)
            classes = mcfit.getclasses()
            results: dict = {k: [] for k in classes.keys()}
            for key, val in classes.items():
                results[key].append(val)
            for mcfit in mcfits:
                classes = mcfit.getclasses()
                for key, val in classes.items():
                    results[key].append(val)
            self._classes = {key: self._getstats(val) for key, val in results.items()}

        return self._classes

    def plot(self, **keywords):
        """
        Plot the MC sampled fit and breakdown components.

        """
        from astropy.nddata import StdDevUncertainty  # type: ignore

        datalabel = keywords.get("datalabel", "obs")

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
                zorder=0,
            )

        else:
            ax.scatter(x, obs.flux.value, color="k", s=5, label=datalabel, zorder=0)

        ax.plot(
            x, fit["mean"], color="tab:purple", linewidth=1.5, label="fit", zorder=100
        )
        ax.fill_between(
            x,
            fit["mean"] - fit["std"],
            fit["mean"] + fit["std"],
            color="tab:purple",
            alpha=0.3,
            zorder=99,
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
                    zorder=2,
                )
                ax.fill_between(
                    x,
                    components["anion"]["mean"] - components["anion"]["std"],
                    components["anion"]["mean"] + components["anion"]["std"],
                    color="tab:red",
                    alpha=0.3,
                    zorder=1,
                )

            if isinstance(components["neutral"]["mean"], np.ndarray):
                ax.plot(
                    x,
                    components["neutral"]["mean"],
                    color="tab:green",
                    linewidth=1.2,
                    label="neutral",
                    zorder=4,
                )
                ax.fill_between(
                    x,
                    components["neutral"]["mean"] - components["neutral"]["std"],
                    components["neutral"]["mean"] + components["neutral"]["std"],
                    color="tab:green",
                    alpha=0.3,
                    zorder=3,
                )

            if isinstance(components["cation"]["mean"], np.ndarray):
                ax.plot(
                    x,
                    components["cation"]["mean"],
                    color="tab:blue",
                    linewidth=1.2,
                    label="cation",
                    zorder=6,
                )
                ax.fill_between(
                    x,
                    components["cation"]["mean"] - components["cation"]["std"],
                    components["cation"]["mean"] + components["cation"]["std"],
                    color="tab:blue",
                    alpha=0.3,
                    zorder=5,
                )

        elif keywords.get("size"):
            ptype = "size"
            if isinstance(components["small"]["mean"], np.ndarray):
                ax.plot(
                    x,
                    components["small"]["mean"],
                    color="tab:red",
                    label="small",
                    zorder=2,
                )
                ax.fill_between(
                    x,
                    components["small"]["mean"] - components["small"]["std"],
                    components["small"]["mean"] + components["small"]["std"],
                    color="tab:red",
                    alpha=0.3,
                    zorder=1,
                )

            if isinstance(components["large"]["mean"], np.ndarray):
                ax.plot(
                    x,
                    components["large"]["mean"],
                    color="tab:green",
                    label="large",
                    zorder=4,
                )
                ax.fill_between(
                    x,
                    components["large"]["mean"] - components["large"]["std"],
                    components["large"]["mean"] + components["large"]["std"],
                    color="tab:green",
                    alpha=0.3,
                    zorder=3,
                )

        elif keywords.get("composition"):
            ptype = "composition"
            if isinstance(components["pure"]["mean"], np.ndarray):
                ax.plot(
                    x,
                    components["pure"]["mean"],
                    color="tab:red",
                    label="pure",
                    zorder=2,
                )
                ax.fill_between(
                    x,
                    components["pure"]["mean"] - components["pure"]["std"],
                    components["pure"]["mean"] + components["pure"]["std"],
                    color="tab:red",
                    alpha=0.3,
                    zorder=1,
                )

            if isinstance(components["nitrogen"]["mean"], np.ndarray):
                ax.plot(
                    x,
                    components["nitrogen"]["mean"],
                    color="tab:green",
                    label="nitrogen",
                    zorder=4,
                )
                ax.fill_between(
                    x,
                    components["nitrogen"]["mean"] - components["nitrogen"]["std"],
                    components["nitrogen"]["mean"] + components["nitrogen"]["std"],
                    color="tab:green",
                    alpha=0.3,
                    zorder=3,
                )

        else:
            ptype = "fitted"

        ax.axhline(0, linestyle="--", color="gray", zorder=0)
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
        if self._error is None:
            mcfits = iter(self.mcfits)
            mcfit = next(mcfits)
            error = mcfit.geterror()
            results: dict = {k: [] for k in error.keys()}
            for key, val in error.items():
                results[key].append(val)
            for mcfit in mcfits:
                error = mcfit.geterror()
                for key, val in error.items():
                    results[key].append(val)
            self._error = {key: self._getstats(val) for key, val in results.items()}

        return self._error

    def getobservation(self) -> Spectrum1D:
        """
        Retrieves the observation.

        """
        return self.observation

    def getresidual(self) -> dict:
        """
        Retrieves the statistics spectra for the residuals of the MC fits.
        """
        if self._residual is None:
            self._residual = self._getstats(
                [mcfit.getresidual() for mcfit in self.mcfits]
            )

        return self._residual
