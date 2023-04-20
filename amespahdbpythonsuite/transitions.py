#!/usr/bin/env python3

from __future__ import annotations

import copy
import multiprocessing
import time
from datetime import timedelta
from functools import partial
from typing import TYPE_CHECKING, Optional

import astropy.units as u  # type: ignore
import numpy as np
from scipy import integrate, optimize  # type: ignore

from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite.data import Data

if TYPE_CHECKING:
    from amespahdbpythonsuite.spectrum import Spectrum


message = AmesPAHdb.message

energy: float
frequency: float
frequencies: np.ndarray
intensities: np.ndarray


class Transitions(Data):
    """
    AmesPAHdbPythonSuite transitions class.
    Contains methods to create and apply an emission model.

    """

    def __init__(self, d: Optional[dict] = None, **keywords) -> None:
        """
        Initialize transitions class.

        """
        super().__init__(d, **keywords)
        self.__set(d, **keywords)

    def set(self, d: Optional[dict] = None, **keywords) -> None:
        """
        Calls class: :class:`amespahdbpythonsuite.data.Data.set` to parse keywords.
        Checks flavor of the database, i.e., theoretical or experimental

        """
        Data.set(self, d, **keywords)
        self.__set(d, **keywords)

    def __set(self, d: Optional[dict] = None, **keywords) -> None:
        """
        Populate data dictionary helper.

        """
        self._shift = keywords.get("shift", 0.0)

        if isinstance(d, dict):
            if d.get("type", "") == self.__class__.__name__:
                if "shift" not in keywords:
                    self._shift = d["shift"]

    def get(self) -> dict:
        """
        Calls class: :class:`amespahdbpythonsuite.data.Data.get` to get keywords.

        """
        d = Data.get(self)
        d["type"] = self.__class__.__name__
        d["shift"] = self._shift
        return copy.deepcopy(d)

    def __repr__(self) -> str:
        """
        Class representation.

        """
        return f"{self.__class__.__name__}(" f"{self.uids=},shift={self._shift})"

    def __str__(self) -> str:
        """
        A description of the instance.

        """

        return f"AmesPAHdbPythonSuite Transitions instance.\n" f"{self.uids=}"

    def write(self, filename: str = "") -> None:
        """
        Write the transitions to file as an IPAC-table.

        """
        import datetime
        import sys

        from astropy.io import ascii  # type: ignore
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
                np.array([t["frequency"] for v in self.data.values() for t in v])
                * self.units["abscissa"]["unit"],
                np.array([t["intensity"] for v in self.data.values() for t in v])
                * self.units["ordinate"]["unit"],
            ],
            names=["UID", "FREQUENCY", "INTENSITY"],
            meta={"comments": hdr},
        )

        if self.database == "theoretical":
            tbl.add_columns(
                [
                    np.array([t["scale"] for v in self.data.values() for t in v]),
                    [t["symmetry"] for v in self.data.values() for t in v],
                ],
                names=["SCALE", "SYMMETRY"],
            )

        ascii.write(tbl, filename, format="ipac", overwrite=True)

        message(f"WRITTEN: {filename}")

    def shift(self, shift: float) -> None:
        """
        Shifts transitions frequency by provided value.

        """
        self._shift += shift

        for key in self.data:
            for d in self.data[key]:
                d["frequency"] += shift
        message(f"TOTAL SHIFT: {self._shift} /cm")

    def fixedtemperature(self, t: float) -> None:
        """
        Applies the Fixed Temperature emission model.

        :param t: Excitation temperature in Kelvin.
        :type t: float

        """
        if self.model:
            if self.model["type"] != "zerokelvin_m":
                message(
                    f'AN EMISSION MODEL HAS ALREADY BEEN APPLIED: {self.model["type"]}'
                )
                return

        message("APPLYING FIXED TEMPERATURE EMISSION MODEL")

        self.model = {
            "type": "fixedtemperature_m",
            "energy": 0.0,
            "temperatures": [t],
            "description": "",
        }

        self.units["ordinate"] = {
            "unit": u.erg / u.second / u.def_unit("molecule", doc="Molecule"),
            "label": "integrated spectral radiance",
        }

        for uid in self.uids:
            f = np.array([d["frequency"] for d in self.data[uid]])

            intensity = (
                2.4853427121856266e-23
                * f**3
                / (np.exp(1.4387751297850830401 * f / t) - 1.0)
            )
            for d, i in zip(self.data[uid], intensity):
                d["intensity"] *= i

    def calculatedtemperature(self, e: float, **keywords) -> None:
        """
        Applies the Calculated Temperature emission model.

        :param e: Excitation energy in erg.
        :type e: float

        """
        if self.database != "theoretical":
            message("THEORETICAL DATABASE REQUIRED FOR EMISSION MODEL")
            return

        if self.model:
            if self.model["type"] != "zerokelvin_m":
                message(
                    f'AN EMISSION MODEL HAS ALREADY BEEN APPLIED: {self.model["type"]}'
                )
                return

        message("APPLYING CALCULATED TEMPERATURE EMISSION MODEL")

        global energy

        energy = e

        self.model = {
            "type": "calculatedtemperature_m",
            "energy": e,
            "temperatures": [],
            "description": "",
        }

        self.units["ordinate"] = {
            "unit": u.erg / u.second / u.def_unit("molecule", doc="Molecule"),
            "label": "integrated spectral radiance",
        }

        print(57 * "=")

        i = 0

        nuids = len(self.uids)

        for uid in self.uids:
            # Start timer.
            tstart = time.perf_counter()

            print("SPECIES                          : %d/%d" % (i + 1, nuids))
            print("UID                              : %d" % uid)
            print(
                "MEAN ABSORBED ENERGY             : %f +/- %f eV"
                % (e / 1.6021765e-12, 0.0)
            )

            global frequencies
            frequencies = np.array([d["frequency"] for d in self.data[uid]])

            Tmax = optimize.brentq(self.attainedtemperature, 2.73, 5000.0)

            print("MAXIMUM ATTAINED TEMPERATURE     : %f Kelvin" % Tmax)

            self.model["temperatures"].append({"uid": uid, "temperature": Tmax})

            for d in self.data[uid]:
                if d["intensity"] > 0:
                    d["intensity"] *= (
                        2.4853427121856266e-23
                        * d["frequency"] ** 3
                        / (np.exp(1.4387751297850830401 * d["frequency"] / Tmax) - 1.0)
                    )

            # Stop timer and calculate elapsed time.
            elapsed = timedelta(seconds=(time.perf_counter() - tstart))
            print(f"Elapsed time: {elapsed}")

            i += 1

        print(57 * "=")

    def cascade(self, e: float, **keywords) -> None:
        """
        Applies the Cascade emission model.

        :param e: Excitation energy in erg.
        :type: float

        """
        if self.database != "theoretical":
            message("THEORETICAL DATABASE REQUIRED FOR EMISSION MODEL")
            return

        if self.model:
            if self.model["type"] != "zerokelvin_m":
                message(
                    f'AN EMISSION MODEL HAS ALREADY BEEN APPLIED: {self.model["type"]}'
                )
                return

        message("APPLYING CASCADE EMISSION MODEL")

        tstart = time.perf_counter()

        global energy

        energy = e

        self.model = {
            "type": "cascade_m",
            "energy": e,
            "temperatures": [],
            "description": "",
        }

        self.units["ordinate"] = {"unit": u.erg, "label": "integrated radiant energy"}

        print(57 * "=")

        if keywords.get("multiprocessing", False):
            cascade_em_model = partial(Transitions._cascade_em_model, e)
            ncores = keywords.get("ncores", multiprocessing.cpu_count() - 1)
            message(f"USING MULTIPROCESSING WITH {ncores} CORES")
            pool = multiprocessing.Pool(processes=ncores)
            data, Tmax = zip(*pool.map(cascade_em_model, self.data.values()))
            pool.close()
            pool.join()

            # Re-assign self.data.
            for uid, d in zip(self.data, data):
                self.data[uid] = d

            self.model["temperatures"] = Tmax

        else:
            i = 0
            nuids = len(self.uids)
            for uid in self.uids:
                print("SPECIES                          : %d/%d" % (i + 1, nuids))
                print("UID                              : %d" % uid)
                print(
                    "MEAN ABSORBED ENERGY             : %f +/- %f eV"
                    % (e / 1.6021765e-12, 0.0)
                )

                global frequencies
                frequencies = np.array([d["frequency"] for d in self.data[uid]])

                global intensities
                intensities = np.array([d["intensity"] for d in self.data[uid]])

                Tmax = optimize.brentq(self.attainedtemperature, 2.73, 5000.0)

                print("MAXIMUM ATTAINED TEMPERATURE     : %f Kelvin" % Tmax)

                self.model["temperatures"].append({"uid": uid, "temperature": Tmax})

                for d in self.data[uid]:
                    if d["intensity"] > 0:
                        global frequency
                        frequency = d["frequency"]
                        d["intensity"] *= (
                            d["frequency"] ** 3
                            * integrate.quad(self.featurestrength, 2.73, Tmax)[0]
                        )

                i += 1

            print(57 * "=")

        elapsed = timedelta(seconds=(time.perf_counter() - tstart))
        print(f"Elapsed time: {elapsed}\n")

    def convolve(self, **keywords) -> Spectrum:
        """
        Convolve transitions with a line profile.
        Calls class:
        :class:`amespahdbpythonsuite.spectrum.Spectrum` to retrieve the
        respective instance.

        """

        fwhm = keywords.get("fwhm", 15.0)

        if keywords.get("gaussian", False):
            width = 0.5 * fwhm / np.sqrt(2.0 * np.log(2.0))
            clip = 3.0
            profile = "Gaussian"
            message("USING GAUSSIAN LINE PROFILES")
        elif keywords.get("drude", False):
            width = 1.0 / fwhm
            clip = 11.0
            profile = "Drude"
            message("USING DRUDE LINE PROFILES")
        else:
            width = 0.5 * fwhm
            clip = 22.0
            profile = "Lorentzian"
            message("USING LORENTZIAN LINE PROFILES")

        if "grid" in keywords:
            x = np.asarray(keywords["grid"])
            xmin = min(x)
            xmax = max(x)
            npoints = len(x)
        else:
            if "xrange" in keywords:
                xmin = min(keywords["xrange"])
                xmax = max(keywords["xrange"])
            else:
                xmin = 1.0
                xmax = 4000.0

            npoints = keywords.get("npoints", 400)
            x = np.arange(xmin, xmax, (xmax - xmin) / npoints)

        message(f"GRID: (XMIN,XMAX)=({xmin:.3f}, {xmax:.3f}); {npoints} POINTS")
        message(f"FWHM: {fwhm} /cm")

        d = dict()

        if keywords.get("multiprocessing", False) and len(self.data) > (
            multiprocessing.cpu_count() - 1
        ):
            get_intensities = partial(
                Transitions._get_intensities,
                npoints,
                xmin,
                xmax,
                clip,
                width,
                x,
                keywords.get("gaussian", False),
                keywords.get("drude", False),
            )
            ncores = keywords.get("ncores", multiprocessing.cpu_count() - 1)
            message(f"USING MULTIPROCESSING WITH {ncores} CORES")
            pool = multiprocessing.Pool(processes=ncores)
            intensities = pool.map(get_intensities, self.data.values())

            pool.close()
            pool.join()

            for uid, i in zip(self.data, intensities):
                d[uid] = i
        else:
            for uid in self.uids:
                s = np.zeros(npoints)
                f = [
                    v
                    for v in self.data[uid]
                    if (v["frequency"] >= xmin - clip * width)
                    and (v["frequency"] <= xmax + clip * width)
                ]
                for t in f:
                    if t["intensity"] > 0:
                        s += t["intensity"] * Transitions._lineprofile(
                            x,
                            t["frequency"],
                            width,
                            gaussian=keywords.get("gaussian", False),
                            drude=keywords.get("drude", False),
                        )
                d[uid] = s

        if self.model["type"] == "zerokelvin_m":
            self.units["ordinate"] = {
                "unit": (u.km * u.cm / u.mol).cgs,
                "label": "cross-section",
            }
        elif (
            self.model["type"] == "fixedtemperature_m"
            or self.model["type"] == "calculatedtemperature_m"
        ):
            self.units["ordinate"] = {
                "unit": u.erg * u.cm / u.mol,
                "label": "radiant energy",
            }
        elif self.model["type"] == "cascade_m":
            self.units["ordinate"] = {
                "unit": u.erg * u.cm / u.mol,
                "label": "radiant energy",
            }

        from amespahdbpythonsuite.spectrum import Spectrum

        return Spectrum(
            database=self.database,
            version=self.version,
            data=d,
            pahdb=self.pahdb,
            uids=self.uids,
            model=self.model,
            units={
                "abscissa": self.units["abscissa"],
                "ordinate": self.units["ordinate"],
            },
            shift=self._shift,
            grid=x,
            profile=profile,
            fwhm=fwhm,
        )

    def plot(self, **keywords) -> None:
        """
        Plot the transitions absorption spectrum.

        """
        import matplotlib.cm as cm  # type: ignore
        import matplotlib.pyplot as plt  # type: ignore

        _, ax = plt.subplots()
        ax.minorticks_on()
        ax.tick_params(which="major", right="on", top="on", direction="in", axis="both")
        colors = cm.rainbow(np.linspace(0, 1, len(self.uids)))
        for uid, col in zip(self.uids, colors):
            f = [v for v in self.data[uid]]
            x = [d["frequency"] for d in f]
            y = [d["intensity"] for d in f]
            ax.bar(x, y, 20, color=col, alpha=0.5)

        plt.gca().invert_xaxis()
        plt.xlabel(
            self.units["abscissa"]["label"]
            + " ["
            + self.units["abscissa"]["unit"].to_string("latex_inline")
            + "]"
        )
        plt.ylabel(
            self.units["ordinate"]["label"]
            + " ["
            + self.units["ordinate"]["unit"].to_string("latex_inline")
            + "]"
        )

        basename = keywords.get("save")
        if basename:
            if not isinstance(basename, str):
                basename = "transitions"
            plt.savefig(f"{basename}.pdf")
        elif keywords.get("show", False):
            plt.show()

    @staticmethod
    def featurestrength(T: float) -> float:
        """
        Calculate a feature's strength covolved with a blackbody.

        :param T: Excitation temperature in Kelvin.
        :type T: float

        """
        global frequency
        global frequencies
        global intensities

        val1 = 1.4387751297850830401 * frequency / T
        if val1 > np.log(np.finfo(float).max):
            return 0.0

        val2 = 1.4387751297850830401 * frequencies / T

        valid = np.where((val2 < np.log(np.finfo(float).max)))

        return (Transitions.heatcapacity(T) / np.expm1(val1)) * (
            1.0
            / np.sum(
                intensities[valid] * (frequencies[valid]) ** 3 / np.expm1(val2[valid])
            )
        )

    @staticmethod
    def attainedtemperature(T: float) -> float:
        """
        Calculate a PAH's temperature after absorbing a given amount of energy.

        :param T: Excitation temperature in Kelvin.
        :type T: float

        """
        global energy

        return integrate.quad(Transitions.heatcapacity, 2.73, T)[0] - energy

    @staticmethod
    def heatcapacity(T: float) -> float:
        """
        Calculate heat capacity.

        :param T: Excitation temperature in Kelvin.
        :type T: float

        """

        global frequencies

        val = 1.4387751297850830401 * frequencies / T

        return 1.3806505e-16 * np.sum(np.exp(-val) * (val / (1.0 - np.exp(-val))) ** 2)

    @staticmethod
    def _lineprofile(x: np.ndarray, x0: float, width: float, **keywords) -> np.ndarray:
        """
        Calculate Gaussian, Drude, or Lorentzian line profiles.

        :param x: Grid array.
        :type x: numpy.ndarray
        :param x0:  Central frequency
        :type x0: float
        :param width: Width of the line profile.
        :type width: float

        """
        if keywords.get("gaussian", False):
            return (1.0 / (width * np.sqrt(2.0 * np.pi))) * np.exp(
                -((x - x0) ** 2) / (2.0 * width**2)
            )
        elif keywords.get("drude", False):
            return (
                (2.0 / (np.pi * x0 * width))
                * width**2
                / ((x / x0 - x0 / x) ** 2 + width**2)
            )
        else:
            return (width / np.pi) / ((x - x0) ** 2 + width**2)

    @staticmethod
    def _get_intensities(
        npoints: int,
        xmin: float,
        xmax: float,
        clip: float,
        width: float,
        x: np.ndarray,
        gaussian: str,
        drude: str,
        data: list,
    ) -> np.ndarray:
        """
        A partial method of :meth:`amespahdbpythonsuite.transitions.convolve`
        used when multiprocessing is required.

        :param npoints: Number of grid points.
        :type npoints: int
        :param xmin: Minimum value of grid.
        :type xmin: float
        :param xmax:  Maximum value of grid.
        :type xmax: float
        :param clip: Value to clip and define the frequency range of the profile
            calculation.
        :type clip: float
        :param width: Width of the line profile.
        :type width: float
        :param x: Grid array
        :type x: numpy.ndarray
        :param gaussian: String to indicate Gaussian profile
        :type gaussian: str
        param drude: String to indicate Drude profile
        type drude: str
        :param data: transitions
        :type data: list

        return s : Transitions
        rtype: numpy.ndarray

        """
        s = np.zeros(npoints)
        f = [
            v
            for v in data
            if v["frequency"] >= xmin - clip * width
            and v["frequency"] <= xmax + clip * width
        ]
        for t in f:
            if t["intensity"] > 0:
                s += t["intensity"] * Transitions._lineprofile(
                    x, t["frequency"], width, gaussian=gaussian, drude=drude
                )

        return s

    @staticmethod
    def _cascade_em_model(e: float, data: list) -> tuple:
        """
        A partial method of :meth:`amespahdbpythonsuite.transitions.cascade`
        used when multiprocessing is required.

        :param e: energy.
        :type e: int

        :param data: transitions.
        :type data: list

        return:  Tupple of transitions and Tmax.
        :rtype: tupple


        """
        global energy
        energy = e

        global frequencies
        frequencies = np.array([d["frequency"] for d in data])

        global intensities
        intensities = np.array([d["intensity"] for d in data])

        Tmax = optimize.brentq(Transitions.attainedtemperature, 2.73, 5000.0)

        for d in data:
            if d["intensity"] > 0:
                global frequency
                frequency = d["frequency"]
                d["intensity"] *= (
                    d["frequency"] ** 3
                    * integrate.quad(Transitions.featurestrength, 2.73, Tmax)[0]
                )

        return data, Tmax
