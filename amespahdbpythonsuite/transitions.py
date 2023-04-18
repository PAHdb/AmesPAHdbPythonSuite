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

    def meanenergyfunc(**keywords):
        # if keywords.get('StellarModel'): 
            # mean_en = 1.9864456023253396e-16 *\
                # integrate.quad(StarModel.frequency, AbsorptionCrossSection__AmesPAHdbIDLSuite(StarModel.frequency) *\
                                    #    StarModel.intensity) /\
                                        # INT_TABULATED(StarModel.frequency, AbsorptionCrosssection__AmesPAHdbIDLSuite(StarModel.frequency) *\
                                                                    # StarModel.intensity / StarModel.frequency)

        # if ketwords.get(ISRF): RETURN,1.9864456023253396D-16 * QROMB('ISRFFunc__AmesPAHdbIDLSuite', 2.5D3, 1.1D5, K=7, EPS=1D-6) / QROMB('ISRFNumberOfPhotonsFunc__AmesPAHdbIDLSuite', 2.5D3, 1.1D5, K=7, EPS=1D-6)

        # RETURN,1.9864456023253396D-16 * QROMB('PlanckFunc__AmesPAHdbIDLSuite', 2.5D3, 1.1D5, K=7, EPS=1D-6) / QROMB('PlanckNumberOfPhotonsFunc__AmesPAHdbIDLSuite', 2.5D3, 1.1D5, K=7, EPS=1D-6)

        return

    def isrf(f: np.ndarray) -> np.ndarray:
        """
        Callback function to calculate the absorption cross-section times
        the interstellar radiation field per Mathis et al. 1983, A&A, 128:212.

        :Params f: frequencies in wavenumber

        :Returns: float array
        """
        # IF f GT 1.1D5 THEN RETURN,0D

        # IF f GT 1D4 / 0.110D THEN RETURN,AbsorptionCrosssection__AmesPAHdbIDLSuite(f) * 1.202D+23 / f^5.4172D

        # IF f GE 1D4 / 0.134D THEN RETURN,AbsorptionCrosssection__AmesPAHdbIDLSuite(f) * 1.366D6 / f^2D

        # IF f GE 1D4 / 0.246D THEN RETURN,AbsorptionCrosssection__AmesPAHdbIDLSuite(f) * 1.019D-2 / f^0.3322D

        # T = [7500D, 4000D, 3000D, 2.73D] & W = [1D-14, 1.65D-13, 4D-13, 1D]

        # RETURN,AbsorptionCrosssection__AmesPAHdbIDLSuite(f) * f^3 * TOTAL(W / (EXP(1.4387751297850830401D * f / T) - 1D))

        return

    def absorptioncrosssection(f: np.ndarray) -> np.ndarray:
        """
        Callback function to Calculate the absorption cross-section
        multiplied by Planck's function.

        :Params f: frequencies in wavenumber

        :Returns: float array
        """
        # ! Need callback for nc, charge, etc. !
        nc = 50
        charge = 1
        # ! Need callback for nc, charge, etc. !

        wave = 1e4 / f

        A_c = [7.97e-17, 1.23e-17, 20e-21, 14e-21, 80e-24, 84e-24, 46e-24, -322e-24]
        W_c = [0.195, 0.217, 0.0805, 0.20, 0.0370, 0.0450, 0.0150, 0.135]
        C_c = [0.0722, 0.2175, 1.05, 1.23, 1.66, 1.745, 1.885, 1.90]

        A = np.transpose(np.resize(A_c, (len(f), len(A_c))))
        W = np.transpose(np.resize(W_c, (len(f), len(W_c))))
        C = np.transpose(np.resize(C_c, (len(f), len(C_c))))

        # Cutoff wavelength from Salama et al. (1996), over wavelength
        # (Eq. (4) in Mattioda et al. (2005))
        y = 1.0 / (0.889 + (2.282 / np.sqrt(0.4 * nc))) / wave

        wave_r2 = np.resize(wave, (2, (len(f))))

        crosssection = ((1.0 / np.pi) * np.arctan((1e3 * (y - 1.0)**3) / y) + 0.5) *\
            (3458e-20 * 10.0**(-3.431 * wave) + (2.0 / np.pi) *\
             np.sum(W[:2] * C[:2] * A[:2] / (((wave_r2 / C[:2]) - (C[:2] / wave_r2))**2 + W[:2]**2), axis=0))

        if charge != 0:

            wave_r6 = np.resize(wave, (6, len(f)))

            crosssection = crosssection + np.exp(-1e-1 / wave**2) * 1.5e-19 * 10**(-wave) + np.sqrt(2.0 / np.pi) *\
                np.sum(A[2:] * np.exp(-2.0 * (wave_r6 - C[2:])**2 / W[2:]**2) / W[2:], axis=0)

        return crosssection

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

        global energy

        energy = e

        tstart = time.perf_counter()

        TStar = 0.

        if keywords.get('Star') and keywords.get('StellarModel'):

            message('STELLAR MODEL SELECTED: USING FIRST PARAMETER AS MODEL')

            TStar = (4 * np.pi * integrate.simpson(e.frequency, e.intensity) / 5.67040e-5)**(0.25)

            select = np.where((e.frequency >= 2.5e3) & (e.frequency <= 1.1e5))
            nselect = e.frequency[select]

            if len(nselect) == 0:
                message('STELLAR MODEL HAS NO DATA BETWEEN 2.5E3-1.1E5 /cm')
                self.state = 0

            message('REBINNING STELLAR MODEL: 100 POINTS')
            StarModel = [{'frequency': 0., 'intensity': 0}] * 100
            sm_frequency = energy['frequency'][select]
            sm_intensity = energy['intensity'][select]
            message(f'CALCULATED EFFECTIVE TEMPERATURE: {TStar} Kelvin')

        elif keywords.get('Star'):
            Tstar = energy
            message(f'BLACKBODY TEMPERATURE: {TStar} Kelvin')

        if keywords.get('ISRF') and not keywords.get('Convolved'):
            message('ISRF SELECTED: IGNORING FIRST PARAMETER')

        if keywords.get('ISRF') or keywords.get('Star') and keywords.get('Convolved'):
            message('CONVOLVING WITH ENTIRE RADIATION FIELD')

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
