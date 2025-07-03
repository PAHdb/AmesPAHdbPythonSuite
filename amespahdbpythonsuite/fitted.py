#!/usr/bin/env python3

import builtins
import operator
import os
from typing import Literal, Optional, Union

import numpy as np
from scipy import integrate  # type: ignore
from specutils import Spectrum1D  # type: ignore

import astropy.units as u  # type: ignore

from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite.spectrum import Spectrum

message = AmesPAHdb.message


class Fitted(Spectrum):
    """
    AmesPAHdbPythonSuite fitted class

    """

    def __init__(self, d: Optional[dict] = None, **keywords) -> None:
        """
        Initialize fitted class.

        """
        super().__init__(d, **keywords)
        self.__set(d, **keywords)

    def plot(self, **keywords) -> None:
        """
        Plotting method for the fitted spectrum and breakdown components.

        """

        import matplotlib as mpl  # type: ignore
        import matplotlib.gridspec as gs  # type: ignore
        import matplotlib.pyplot as plt  # type: ignore

        if not self.observation:
            message("OBSERVATIONS ARE MISSING")
            return

        if keywords.get("datalabel", False):
            datalabel = keywords["datalabel"]
        else:
            datalabel = "obs"

        axis = []
        if keywords.get("sizedistribution", False):
            fig, ax = plt.subplots()
            axis.append(ax)
            h, edges = self.getsizedistribution()
            h = 100.0 * h / np.sum(h)
            axis[0].bar(
                edges[:-1],
                h,
                align="edge",
                edgecolor="black",
                width=(np.roll(edges, -1) - edges)[:-1],
            )
            axis[0].set_xlabel(r"n$_{\mathregular{carbon}}$")
            axis[0].set_ylabel("frequency [%]")
        else:
            if keywords.get("residual", False):
                fig = plt.figure()
                spec = gs.GridSpec(2, 1, height_ratios=[3, 1])
                axis.append(plt.subplot(spec[0]))
                axis.append(plt.subplot(spec[1], sharex=axis[0]))
                fig.subplots_adjust(hspace=0)
                axis[0].tick_params(axis="x", which="both", labelbottom="off")
            elif (
                keywords.get("charge", False)
                or keywords.get("size", False)
                or keywords.get("composition", False)
            ):
                fig, ax = plt.subplots()
                axis.append(ax)
            elif not keywords.get("sizedistribution", False):
                fig = plt.figure()
                spec = gs.GridSpec(1, 2, width_ratios=[2, 3])
                axis.append(plt.subplot(spec[0]))
                axis.append(plt.subplot(spec[1]))
                fig.subplots_adjust(wspace=0.25)
                axis[0].tick_params(
                    axis="x", which="both", bottom="off", top="off", labelbottom="off"
                )
                axis[0].tick_params(
                    axis="y", which="both", left="off", right="off", labelleft="off"
                )
                axis[0].set_ylim((0, 1))
                axis = list(reversed(axis))

            if keywords.get("wavelength", False):
                x = 1e4 / self.grid
                axis[0].set_xlim((min(x), max(x)))
                xtitle = "Wavelength [$\\mu$m]"
            else:
                x = self.grid
                axis[0].set_xlim((max(x), min(x)))
                xtitle = (
                    self.units["abscissa"]["label"]
                    + " ["
                    + self.units["abscissa"]["unit"].to_string("latex_inline")
                    + "]"
                )

            if "sigma" in keywords:
                axis[0].errorbar(
                    x,
                    self.observation.flux.value,
                    yerr=keywords["sigma"],
                    fmt="o",
                    mfc="white",
                    color="k",
                    ecolor="k",
                    markersize=2,
                    elinewidth=0.2,
                    capsize=0.8,
                    label=datalabel,
                )
            else:
                axis[0].scatter(
                    x, self.observation.flux.value, color="k", s=1, label=datalabel
                )

            if "title" in keywords:
                axis[0].set_title(keywords["title"])

            axis[0].minorticks_on()
            axis[0].tick_params(
                which="major", right="on", top="on", direction="in", length=5
            )
            axis[0].tick_params(
                which="minor", right="on", top="on", direction="in", length=3
            )

            if (
                not keywords.get("residual", False)
                and not keywords.get("charge", False)
                and not keywords.get("size", False)
                and not keywords.get("composition", False)
            ):
                colors = mpl.colormaps["rainbow"](np.linspace(0, 1, len(self.uids)))
                for uid, col in zip(self.uids, colors):
                    axis[0].plot(x, self.data[uid], color=col)

            axis[0].plot(
                x, self.getfit(), color="tab:purple", label="fit", lw=1, zorder=42
            )

            if (
                keywords.get("charge", False)
                or keywords.get("size", False)
                or keywords.get("composition", False)
            ):
                classes = self.getclasses()
                if isinstance(classes, dict):
                    if keywords.get("charge", False):
                        if not isinstance(classes["anion"], int):
                            axis[0].plot(
                                x, classes["anion"], color="tab:red", label="anion"
                            )
                        if not isinstance(classes["neutral"], int):
                            axis[0].plot(
                                x,
                                classes["neutral"],
                                color="tab:green",
                                label="neutral",
                            )
                        if not isinstance(classes["cation"], int):
                            axis[0].plot(
                                x, classes["cation"], color="tab:blue", label="cation"
                            )
                        axis[0].axhline(0, linestyle="--", color="gray")
                        axis[0].legend()
                    elif keywords.get("size", False):
                        if not isinstance(classes["small"], int):
                            axis[0].plot(
                                x, classes["small"], color="tab:blue", label="small"
                            )
                        if not isinstance(classes["medium"], int):
                            axis[0].plot(
                                x, classes["medium"], color="tab:green", label="medium"
                            )
                        if not isinstance(classes["large"], int):
                            axis[0].plot(
                                x, classes["large"], color="tab:red", label="large"
                            )
                        axis[0].axhline(0, linestyle="--", color="gray")
                        axis[0].legend()
                    elif keywords.get("composition", False):
                        if not isinstance(classes["pure"], int):
                            axis[0].plot(
                                x, classes["pure"], color="tab:red", label="pure"
                            )
                        if not isinstance(classes["nitrogen"], int):
                            axis[0].plot(
                                x,
                                classes["nitrogen"],
                                color="tab:green",
                                label="nitrogen",
                            )
                        axis[0].axhline(0, linestyle="--", color="gray")
                        axis[0].legend()
            elif keywords.get("residual", False):
                y = self.getresidual()
                if y is not None:
                    axis[1].plot(x, y)
                axis[1].axhline(0, linestyle="--", color="gray")
                axis[0].legend()
            else:
                axis[1].text(
                    0.05,
                    0.95,
                    ("%s" + 5 * " " + "%s") % ("UID", "WEIGHT"),
                    family="monospace",
                )
                axis[1].xaxis.set_visible(False)
                axis[1].yaxis.set_visible(False)
                ypos = 0.95 - 0.05
                for uid, w, col in zip(self.uids, self.weights.values(), colors):
                    txt = ("%d" + (5 * " ") + "%.2e") % (uid, w)
                    axis[1].text(0.05, ypos, txt, color=col, family="monospace")
                    ypos -= 0.05
                    if ypos <= 0.05:
                        axis[1].text(0.05, ypos, "more...", family="monospace")
                        break

            axis[0].ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))

            if self.observation.flux.unit == u.Unit():
                axis[0].set_ylabel(
                    self.units["ordinate"]["label"]
                    + " ["
                    + self.units["ordinate"]["unit"].to_string("latex_inline")
                    + "]",
                )
            else:
                axis[0].set_ylabel(self.observation.flux.unit.to_string("latex_inline"))
            if keywords.get("residual", False):
                axis[1].set_xlabel(f"{xtitle}")
                axis[1].set_ylabel("residual")
                axis[1].minorticks_on()
                axis[1].tick_params(
                    which="major", right="on", top="on", direction="in", length=5
                )
                axis[1].tick_params(
                    which="minor", right="on", top="on", direction="in", length=3
                )
            else:
                axis[0].set_xlabel(f"{xtitle}")

        if keywords.get("save", False):
            if keywords.get("charge"):
                ptype = "charge"
            elif keywords.get("size"):
                ptype = "size"
            elif keywords.get("composition"):
                ptype = "composition"
            elif keywords.get("residual"):
                ptype = "residual"
            else:
                ptype = "fitted"

            if keywords["output"]:
                if os.path.isdir(keywords["output"]):
                    fig.savefig(
                        f"{keywords['output']}/{ptype}.{keywords['ftype']}",
                        bbox_inches="tight",
                    )
                else:
                    fig.savefig(
                        f"{keywords['output']}_{ptype}.{keywords['ftype']}",
                        bbox_inches="tight",
                    )
            else:
                fig.savefig(f"{ptype}.{keywords['ftype']}", bbox_inches="tight")
        else:
            plt.show()
        plt.close(fig)

    def set(self, d: Optional[dict] = None, **keywords) -> None:
        """
        Calls class: :class:`amespahdbpythonsuite.spectrum.Spectrum.set` to parse keywords.

        """
        Spectrum.set(self, d, **keywords)
        self.__set(d, **keywords)

    def __set(self, d: Optional[dict] = None, **keywords) -> None:
        """
        Populate data dictionary helper.

        """
        self.observation = keywords.get("observation", None)
        self.weights = keywords.get("weights", dict())
        self.method = keywords.get("method", "")

        if isinstance(d, dict):
            if d.get("type", "") == self.__class__.__name__:
                if "observation" not in keywords:
                    self.observation = d["observation"]
                if "weights" not in keywords:
                    self.weights = d["weights"]
                if "method" not in keywords:
                    self.method = d["method"]

        self._chisquared: Optional[float] = None
        self._norm: Optional[float] = None
        self._fit: Optional[np.ndarray] = None
        self._classes: Optional[dict] = None
        self._error: Optional[dict] = None
        self._residual: Optional[np.ndarray] = None

    def get(self) -> dict:
        """
        Calls class: :class:`amespahdbpythonsuite.spectrum.Spectrum.get`.
        Assigns class variables from inherited dictionary.

        """
        d = Spectrum.get(self)
        d["type"] = self.__class__.__name__
        d["observation"] = self.observation
        d["weights"] = self.weights
        d["method"] = self.method

        return d

    def __repr__(self) -> str:
        """
        Class representation.

        """
        return f"{self.__class__.__name__}({self.uids=},{self.method=})"

    def __str__(self) -> str:
        """
        A description of the instance.
        """

        return f"AmesPAHdbPythonSuite Fitted instance.\n{self.uids=}"

    def write(self, filename: str = "") -> None:
        """
        Write the fitted spectra to file as an IPAC-table.

        """
        import datetime
        import sys

        from astropy.io import ascii  # type: ignore
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
            "SPECIES": str(len(self.data)),
        }

        for key, value in kv.items():
            if not value.isnumeric():
                hdr.append(f"{key:8} = '{value}'")
            else:
                hdr.append(f"{key:8} = {value}")

        tbl = Table(
            [
                [uid for uid, v in self.data.items() for _ in v],
                np.array([f for _ in self.data.values() for f in self.grid])
                * self.units["abscissa"]["unit"],
                np.array([t for v in self.data.values() for t in v])
                * self.units["ordinate"]["unit"],
                np.array(
                    [self.weights[uid] for uid, v in self.data.items() for _ in v]
                ),
            ],
            names=["UID", "FREQUENCY", "INTENSITY", "WEIGHT"],
            meta={"comments": hdr},
        )

        ascii.write(tbl, filename, format="ipac", overwrite=True)

        message(f"WRITTEN: {filename}")

    def getmethod(self) -> str:
        """
        Retrieves the method used for the fit.

        """
        return self.method

    def getchisquared(self) -> Optional[float]:
        """
        Calculates the chi-squared of the fit.

        """
        if self._chisquared is None:
            if self.observation:
                self._chisqured = np.sum(
                    (self.observation.flux.value - self.getfit()) ** 2
                    / self.observation.uncertainty.array
                )

        return self._chisquared

    def getnorm(self) -> Optional[float]:
        """
        Retrieves the Norm of the fit.

        """
        if self._norm is None and self.observation:
            self._norm = np.sum((self.observation.flux.value - self.getfit()) ** 2)

        return self._norm

    def getobservation(self) -> Spectrum1D:
        """
        Retrieves the ordinate values of the observation.

        """
        return self.observation

    def getfit(self) -> Union[np.ndarray, Literal[0]]:
        """
        Retrieves the sum of fit values.

        """
        if self._fit is None:
            self._fit = sum(self.data.values())

        return self._fit

    def getresidual(self) -> Optional[np.ndarray]:
        """
        Retrieves the residual of the fit.
        """
        if self._residual is None and self.observation:
            self._residual = self.observation.flux.value - self.getfit()

        return self._residual

    def getweights(self) -> dict:
        """
        Retrieves the weights of the fitted PAHs.

        """
        return self.weights

    def getsizedistribution(
        self,
        nbins: int = 0,
        min: Optional[float] = None,
        max: Optional[float] = None,
    ) -> tuple:
        """
        Retrieves the size distribution of the fitted PAHs.

        """

        if not self.pahdb:
            return ()

        nc = [self.pahdb["species"][uid]["n_c"] for uid in self.weights]

        if not min:
            min = builtins.min(nc)
        if not max:
            max = builtins.max(nc)
        if nbins == 0:
            nbins = int(np.ceil(2.0 * len(self.uids) ** (1.0 / 3.0)))

        return np.histogram(nc, bins=nbins, weights=list(self.weights.values()))

    def getclasses(self, **keywords) -> Optional[dict]:
        """
        Retrieves the spectra of the different classes of the fitted PAHs.

        """
        if self._classes is None:
            if self.pahdb is None:
                message("VALID DATABASE NEEDED TO GET CLASSES")
                return None

            self._classes = {
                key: self.__classes(val)
                for key, val in self._subclasses(**keywords).items()
            }

            uids = [
                uid
                for uid in self.uids
                if self.pahdb["species"][uid]["n_n"] == 0
                and self.pahdb["species"][uid]["n_o"] == 0
                and self.pahdb["species"][uid]["n_mg"] == 0
                and self.pahdb["species"][uid]["n_si"] == 0
                and self.pahdb["species"][uid]["n_fe"] == 0
            ]
            self._classes["pure"] = sum(
                {key: val for key, val in self.data.items() if key in uids}.values()
            )

        return self._classes

    def __classes(self, s: dict) -> Optional[np.ndarray]:
        """
        Retrieves the intensities of a given subclass.

        """

        if not self.pahdb:
            return None

        if s["subclass"] == "n_c":
            if "operator_1" in s:
                uids = [
                    uid
                    for uid in self.uids
                    if s["operator_1"](
                        self.pahdb["species"][uid][s["subclass"]],
                        s["operand_1"],
                    )
                    and s["operator_2"](
                        self.pahdb["species"][uid][s["subclass"]],
                        s["operand_2"],
                    )
                ]
            else:
                uids = [
                    uid
                    for uid in self.uids
                    if s["operator"](
                        self.pahdb["species"][uid][s["subclass"]], s["operand"]
                    )
                ]
        else:
            uids = [
                uid
                for uid in self.uids
                if s["operator"](
                    self.pahdb["species"][uid][s["subclass"]], s["operand"]
                )
            ]

        if not uids:
            return np.zeros(self.grid.size)

        return sum({key: val for key, val in self.data.items() if key in uids}.values())

    def getbreakdown(self, **keywords) -> Optional[dict]:
        """
        Retrieves the breakdown of the fitted PAHs.

        """
        if not self.pahdb:
            message("VALID DATABASE NEEDED TO GET CLASSES")
            return None

        # Getting fit weights
        fweight = np.array(list(self.weights.values()))

        breakdown = {
            "solo": np.sum([self.pahdb["species"][uid]["n_solo"] for uid in self.uids]),
            "duo": np.sum([self.pahdb["species"][uid]["n_duo"] for uid in self.uids]),
            "trio": np.sum([self.pahdb["species"][uid]["n_trio"] for uid in self.uids]),
            "quartet": np.sum(
                [self.pahdb["species"][uid]["n_quartet"] for uid in self.uids]
            ),
            "quintet": np.sum(
                [self.pahdb["species"][uid]["n_quintet"] for uid in self.uids]
            ),
        }

        total = 1.0

        if keywords.get("flux", False):
            if not self.grid:
                message("GRID IS NOT SET")
                return None

            classes = self.getclasses(**keywords)

            if keywords.get("absolute", False):
                total, err = -integrate.trapezoid(self.getfit(), x=self.grid)

            if classes:
                for key in classes:
                    breakdown[key] = (
                        -integrate.trapezoid(classes[key], x=self.grid) / total
                    )

            return breakdown

        if not self.weights:
            message("WEIGHTS ARE NOT SET")

            return None

        if "absolute" not in keywords:
            total = np.sum(fweight)

        # Set subclasses dictionary.
        subclasses = self._subclasses(**keywords)

        for key in subclasses:
            breakdown[key] = self.__breakdown(subclasses[key]) / total

        # Obtaining pure PAH breakdown.
        uids = [
            uid
            for uid in self.uids
            if self.pahdb["species"][uid]["n_n"] == 0
            and self.pahdb["species"][uid]["n_o"] == 0
            and self.pahdb["species"][uid]["n_mg"] == 0
            and self.pahdb["species"][uid]["n_si"] == 0
            and self.pahdb["species"][uid]["n_fe"] == 0
        ]

        if len(uids) > 0:
            breakdown["pure"] = np.sum([self.weights[uid] for uid in uids]) / total

        # Getting Nc.
        nc = np.array([self.pahdb["species"][uid]["n_c"] for uid in self.uids])
        breakdown["n_c"] = np.sum(nc * fweight) / np.sum(fweight)

        return breakdown

    def __breakdown(self, s: dict) -> float:
        """
        Retrieve the sum of the fitting weights for the fitted PAHs.

        """
        if not self.pahdb:
            return 0.0

        if "operator_1" in s:
            uids = [
                uid
                for uid in self.uids
                if s["operator_1"](
                    self.pahdb["species"][uid][s["subclass"]],
                    s["operand_1"],
                )
                and s["operator_2"](
                    self.pahdb["species"][uid][s["subclass"]],
                    s["operand_2"],
                )
            ]
        else:
            uids = [
                uid
                for uid in self.uids
                if s["operator"](
                    self.pahdb["species"][uid][s["subclass"]], s["operand"]
                )
            ]

        return np.sum([self.weights[uid] for uid in uids]) if len(uids) else 0.0

    def _subclasses(self, **keywords) -> dict[str, dict]:
        """
        Create subclasses dictionary.

        """

        return {
            "anion": {"subclass": "charge", "operator": operator.lt, "operand": 0},
            "neutral": {"subclass": "charge", "operator": operator.eq, "operand": 0},
            "cation": {"subclass": "charge", "operator": operator.gt, "operand": 0},
            "small": {
                "subclass": "n_c",
                "operator": operator.le,
                "operand": keywords.get("small", 50),
            },
            "medium": {
                "subclass": "n_c",
                "operator_1": operator.gt,
                "operator_2": operator.le,
                "operand_1": keywords.get("medium_low", 50),
                "operand_2": keywords.get("medium_high", 70),
            },
            "large": {
                "subclass": "n_c",
                "operator": operator.gt,
                "operand": keywords.get("large", 70),
            },
            "nitrogen": {"subclass": "n_n", "operator": operator.gt, "operand": 0},
        }

    def geterror(self) -> Optional[dict[str, float]]:
        """
        Calculates the PAHdb fitting uncertainty
        as the ratio of the residual over the total spectrum area.

        """
        if self._error is None:
            range = {
                "err": [min(self.grid), max(self.grid)],
                "e127": [754.0, 855.0],
                "e112": [855.0, 1000.0],
                "e77": [1000.0, 1495.0],
                "e62": [1495.0, 1712.0],
                "e33": [2900.0, 3125.0],
            }

            self._error = {k: 0.0 for k in range.keys()}
            if self.observation:
                srt = np.argsort(self.grid)
                x = self.grid[srt]
                y = self.observation.flux.value[srt]
                r = self.getresidual()
                if r is not None:
                    r = np.abs(r)[srt]
                    if isinstance(r, np.ndarray):
                        for key, rng in range.items():
                            sel = np.where((x >= rng[0]) & (x <= rng[1]))[0]

                            if len(sel) == 0:
                                continue

                            self._error[key] = np.divide(
                                integrate.trapezoid(r[sel], x=x[sel]),
                                integrate.trapezoid(y[sel], x=x[sel]),
                            )

        return self._error

    def sort(self, flux: bool = False) -> dict:
        """
        Sort UIDs and weights by their contribution to the fit

        """
        w = (
            self.weights
            if not flux
            else {
                uid: -integrate.trapezoid(flux, x=self.grid)
                for uid, flux in self.data.items()
            }
        )

        self.weights = dict(
            sorted(w.items(), key=lambda item: (item[1], item[0]), reverse=True)
        )

        self.uids = list(self.weights.keys())

        return self.weights
