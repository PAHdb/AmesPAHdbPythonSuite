#!/usr/bin/env python3

from typing import Optional, Union, Literal

import os
import builtins
import operator
import numpy as np

from scipy import integrate  # type: ignore
from specutils import Spectrum1D  # type: ignore

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
        self.atoms: dict = dict()
        self.__set(d, **keywords)

    def plot(self, **keywords) -> None:
        """
        Plotting method for the fitted spectrum and breakdown components.

        """

        import matplotlib.pyplot as plt  # type: ignore
        import matplotlib.gridspec as gs  # type: ignore
        import matplotlib.cm as cm  # type: ignore

        if keywords.get("datalabel", False):
            datalabel = keywords["datalabel"]
        else:
            datalabel = "obs"

        if keywords.get("sizedistribution", False):
            fig, axis = plt.subplots()
            axis = [axis]
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
                axis = plt.subplot(spec[0])
                axis = [axis, plt.subplot(spec[1], sharex=axis)]
                fig.subplots_adjust(hspace=0)
                axis[0].tick_params(axis="x", which="both", labelbottom="off")
            elif (
                keywords.get("charge", False)
                or keywords.get("size", False)
                or keywords.get("composition", False)
            ):
                fig, axis = plt.subplots()
                axis = [axis]
            elif not keywords.get("sizedistribution", False):
                fig = plt.figure()
                spec = gs.GridSpec(1, 2, width_ratios=[2, 3])
                axis = plt.subplot(spec[0])
                axis = [axis, plt.subplot(spec[1])]
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
                    markersize=3,
                    elinewidth=0.2,
                    capsize=0.8,
                    label=datalabel,
                )
            else:
                axis[0].scatter(
                    x, self.observation.flux.value, color="k", s=5, label=datalabel
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

            colors = cm.rainbow(np.linspace(0, 1, len(self.uids)))
            for uid, col in zip(self.uids, colors):
                y = self.data[uid]
                if keywords.get("residual", False):
                    axis[0].plot(x, y, color=col)

            axis[0].plot(x, self.getfit(), color="tab:purple", label="fit")

            if (
                keywords.get("charge", False)
                or keywords.get("size", False)
                or keywords.get("composition", False)
            ):
                classes = self.getclasses()
                if keywords.get("charge", False):
                    if not isinstance(classes["anion"], int):
                        axis[0].plot(x, classes["anion"], "r-", label="anion")
                    if not isinstance(classes["neutral"], int):
                        axis[0].plot(x, classes["neutral"], "g-", label="neutral")
                    if not isinstance(classes["cation"], int):
                        axis[0].plot(x, classes["cation"], "b-", label="cation")
                    axis[0].axhline(0, linestyle="--", color="gray")
                    axis[0].legend()
                elif keywords.get("size", False):
                    if not isinstance(classes["small"], int):
                        axis[0].plot(x, classes["small"], "r-", label="small")
                    if not isinstance(classes["large"], int):
                        axis[0].plot(x, classes["large"], "g-", label="large")
                    axis[0].axhline(0, linestyle="--", color="gray")
                    axis[0].legend()
                elif keywords.get("composition", False):
                    if not isinstance(classes["pure"], int):
                        axis[0].plot(x, classes["pure"], "r-", label="pure")
                    if not isinstance(classes["nitrogen"], int):
                        axis[0].plot(x, classes["nitrogen"], "g-", label="nitrogen")
                    axis[0].axhline(0, linestyle="--", color="gray")
                    axis[0].legend()
            elif keywords.get("residual", False):
                y = self.getresidual()
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

            axis[0].set_ylabel(
                self.units["ordinate"]["label"]
                + " ["
                + self.units["ordinate"]["unit"].to_string("latex_inline")
                + "]",
            )
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
                    fig.savefig(f"{keywords['output']}/{ptype}.{keywords['ftype']}")
                else:
                    fig.savefig(f"{keywords['output']}_{ptype}.{keywords['ftype']}")
            else:
                fig.savefig(f"{ptype}.{keywords['ftype']}")
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
        return f"{self.__class__.__name__}(" f"{self.uids=},{self.method=})"

    def __str__(self) -> str:
        """
        A description of the instance.
        """

        return f"AmesPAHdbPythonSuite Fitted instance.\n" f"{self.uids=}"

    def write(self, filename: str = "") -> None:
        """
        Write the fitted spectra to file as an IPAC-table.

        """
        import sys
        import datetime
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
        if self.observation:
            return np.sum(
                (self.observation.flux.value - self.getfit()) ** 2
                / self.observation.uncertainty.array
            )

        return None

    def getnorm(self) -> float:
        """
        Retrieves the Norm of the fit.

        """
        return np.sum((self.observation.flux.value - self.getfit()) ** 2)

    def getobservation(self) -> Spectrum1D:
        """
        Retrieves the ordinate values of the observation.

        """
        return self.observation

    def getfit(self) -> Union[np.ndarray, Literal[0]]:
        """
        Retrieves the sum of fit values.

        """
        return sum(self.data.values())

    def getresidual(self) -> float:
        """
        Retrieves the residual of the fit.
        """
        return self.observation.flux.value - self.getfit()

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

        if not len(self.atoms):
            self._atoms()

        nc = [self.atoms[uid]["nc"] for uid in self.weights]

        if not min:
            min = builtins.min(nc)
        if not max:
            max = builtins.max(nc)
        if nbins == 0:
            nbins = int(np.ceil(2.0 * len(self.uids) ** (1.0 / 3.0)))

        return np.histogram(nc, bins=nbins, weights=list(self.weights.values()))

    def getclasses(self, **keywords) -> dict:
        """
        Retrieves the spectra of the different classes of the fitted PAHs.

        """
        if not self.pahdb:
            message("VALID DATABASE NEEDED TO GET CLASSES")
            return dict()

        if not len(self.atoms):
            self._atoms()

        # Set subclasses dictionary.
        subclasses = self._subclasses(**keywords)

        classes: dict = dict()

        for key in subclasses:
            classes[key] = self.__classes(subclasses[key])

        puids = [
            uid
            for uid in self.uids
            if self.atoms[uid]["nn"] == 0
            and self.atoms[uid]["no"] == 0
            and self.atoms[uid]["nmg"] == 0
            and self.atoms[uid]["nsi"] == 0
            and self.atoms[uid]["nfe"] == 0
        ]
        classes["pure"] = sum(
            {k: v for k, v in self.data.items() if k in puids}.values()
        )

        return classes

    def __classes(self, s: dict) -> np.ndarray:
        """
        Retrieves the intensities of a given subclass.

        """

        if s["subclass"] == "charge":
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
                if s["operator"](self.atoms[uid][s["subclass"]], s["operand"])
            ]

        if not uids:
            return np.zeros(self.grid.size)

        return sum({k: v for k, v in self.data.items() if k in uids}.values())

    def getbreakdown(self, **keywords) -> Optional[dict]:
        """
        Retrieves the breakdown of the fitted PAHs.

        """
        if not self.pahdb:
            message("VALID DATABASE NEEDED TO GET CLASSES")
            return None

        # Set atom dictionary if it doesn't exist.
        if not len(self.atoms):
            self._atoms()

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
            if self.atoms[uid]["nn"] == 0
            and self.atoms[uid]["no"] == 0
            and self.atoms[uid]["nmg"] == 0
            and self.atoms[uid]["nsi"] == 0
            and self.atoms[uid]["nfe"] == 0
        ]

        if len(uids) > 0:
            breakdown["pure"] = np.sum([self.weights[uid] for uid in uids]) / total

        # Getting Nc.
        nc = np.array([self.atoms[uid]["nc"] for uid in self.uids])
        breakdown["nc"] = np.sum(nc * fweight) / np.sum(fweight)

        return breakdown

    def __breakdown(self, s: dict) -> float:
        """
        Retrieve the sum of the fitting weights for the fitted PAHs.

        """
        if s["subclass"] == "charge":
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
                if s["operator"](self.atoms[uid][s["subclass"]], s["operand"])
            ]

        if len(uids) > 0:
            return np.sum([self.weights[uid] for uid in uids])
        else:
            return 0.0

    def _subclasses(self, **keywords) -> dict:
        """
        Create subclasses dictionary.

        """

        subclasses = {
            "anion": {"subclass": "charge", "operator": operator.lt, "operand": 0},
            "neutral": {"subclass": "charge", "operator": operator.eq, "operand": 0},
            "cation": {"subclass": "charge", "operator": operator.gt, "operand": 0},
            "small": {
                "subclass": "nc",
                "operator": operator.le,
                "operand": keywords.get("small", 50),
            },
            "large": {
                "subclass": "nc",
                "operator": operator.gt,
                "operand": keywords.get("small", 50),
            },
            "nitrogen": {"subclass": "nn", "operator": operator.gt, "operand": 0},
        }

        return subclasses

    def _atoms(self) -> None:
        """
        Create atoms dictionary.

        """

        # Create reference dictionary with atomic numbers for c, h, n, o, mg, si, and fe.
        nelem = {"nc": 6, "nh": 1, "nn": 7, "no": 8, "nmg": 12, "nsi": 14, "nfe": 26}

        # Initialize atoms dictionary.
        self.atoms = {key: {} for key in self.uids}

        for uid in self.uids:
            # Initialize dictionary based on reference dictionary.
            dnelem = dict.fromkeys(nelem)
            # Loop through the keys to determine the number of each atom present in a given uid.
            for key in dnelem.keys():
                dnelem[key] = len(
                    [
                        sub["type"]
                        for sub in self.pahdb["species"][uid]["geometry"]
                        if sub["type"] == nelem[key]
                    ]
                )
            self.atoms[uid] = dnelem

    def geterror(self) -> Optional[dict]:
        """
        Calculates the PAHdb fitting uncertainty
        as the ratio of the residual over the total spectrum area.

        """
        tags = ["err", "e127", "e112", "e77", "e62", "e33"]

        piecewise = dict.fromkeys(tags, None)

        range = dict()
        range["err"] = [min(self.grid), max(self.grid)]
        range["e127"] = [754.0, 855.0]
        range["e112"] = [855.0, 1000.0]
        range["e77"] = [1000.0, 1495.0]
        range["e62"] = [1495.0, 1712.0]
        range["e33"] = [2900.0, 3125.0]

        if self.observation:
            for key in piecewise.keys():
                sel = np.where(
                    np.logical_and(
                        self.grid >= range[key][0], self.grid <= range[key][1]
                    )
                )[0]
                total_area = np.trapz(
                    self.observation.flux.value[sel], x=self.grid[sel]
                )
                if total_area == 0:
                    continue
                fit = self.getfit()
                if isinstance(fit, np.ndarray):
                    resid_area = np.trapz(
                        np.abs(self.observation.flux.value[sel] - fit[sel]),
                        x=self.grid[sel],
                    )
                    piecewise[key] = resid_area / total_area

        return piecewise

    def sort(self, flux: bool = False) -> dict:
        """
        Sort UIDs and weights by their contribution to the fit

        """

        if not flux:
            w = self.weights
        else:
            w = dict()
            for uid, flux in self.data.items():
                w[uid] = -integrate.trapezoid(flux, x=self.grid)

        self.weights = dict(
            sorted(w.items(), key=lambda item: (item[1], item[0]), reverse=True)
        )

        self.uids = list(self.weights.keys())

        return self.weights
