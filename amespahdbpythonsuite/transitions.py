#!/usr/bin/env python3

from __future__ import annotations

import copy
import multiprocessing
import time
from datetime import timedelta
from functools import partial
from typing import TYPE_CHECKING, Optional

from scipy import integrate, optimize, ndimage, special  # type: ignore

import astropy.units as u  # type: ignore
import numpy as np

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

    def fixed_temperature(self, t: float) -> None:
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
                / (np.expm1(1.4387751297850830401 * f / t))
            )
            for d, i in zip(self.data[uid], intensity):
                d["intensity"] *= i

    def calculated_temperature(self, e: float, **keywords) -> None:
        """
        Applies the Calculated Temperature emission model.

        :param e: Excitation energy in erg.
        :type e: float
        :If 'star' keyword is given the input should be temperature in Kelvin instead of energy.

        :keywords:
        :

        """
        if self.database != "theoretical":
            message("THEORETICAL DATABASE REQUIRED FOR EMISSION MODEL")
            return

        # if keywords.get('star') or keywords.get('isrf') and not self.database:
        # message("VALID DATABASE NEEDED FOR USE WITH STAR/ISRF")
        # return

        if self.model:
            if self.model["type"] != "zerokelvin_m":
                message(
                    f'AN EMISSION MODEL HAS ALREADY BEEN APPLIED: {self.model["type"]}'
                )
                return

        message("APPLYING CALCULATED TEMPERATURE EMISSION MODEL")

        global energy
        global Tstar

        energy = e
        Tstar = 0

        if (keywords.get('star')) and (keywords.get('stellar_model')):
            message('STELLAR MODEL SELECTED: USING FIRST PARAMETER AS MODEL')

            Tstar = (
                4
                * np.pi
                * integrate.simpson(energy["frequency"], energy["intensity"])
                / 5.67040e-5
            ) ** 0.25

            select = np.where(
                (energy["frequency"] >= 2.5e3) & (energy["frequency"] <= 1.1e5)
            )
            nselect = energy["frequency"][select]

            if len(nselect) == 0:
                message('STELLAR MODEL HAS NO DATA BETWEEN 2.5E3-1.1E5 /cm')
                self.state = 0

            message('REBINNING STELLAR MODEL: 100 POINTS')
            # ^ Consider giving a user warning when providing large dimension data.

            global star_model
            star_model = [{'frequency': 0.0, 'intensity': 0}] * 100
            # * Using interpolation through ndimage.zoom to get equivalent
            # * of the IDL congrid function.
            star_model['frequency'] = ndimage.zoom(
                energy['frequency'][select], 100 / len(energy['frequency'][select])
            )
            star_model['intensity'] = ndimage.zoom(
                energy['intensity'][select], 100 / len(energy['intensity'][select])
            )
            message(f'CALCULATED EFFECTIVE TEMPERATURE: {Tstar} Kelvin')

        elif keywords.get('star'):
            Tstar = energy

        if keywords.get('isrf') and not keywords.get('convolved'):
            message('ISRF SELECTED: IGNORING FIRST PARAMETER')

        self.model = {
            "type": "calculated_temperature_m",
            "energy": {},
            "temperatures": {},
            "approximate": keywords.get('approximate'),
            "star": keywords.get('star'),
            "isrf": keywords.get('isrf'),
            "stellar_model": keywords.get('stellar_model'),
            "Tstar": Tstar,
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
            # Instatiate model energy and temperature dictionaries.
            self.model['energy'][uid] = {}
            self.model["temperatures"][uid] = {}

            print("SPECIES                          : %d/%d" % (i + 1, nuids))
            print("UID                              : %d" % uid)

            if (
                keywords.get('approximate')
                or keywords.get('star')
                or keywords.get('isrf')
            ):
                global charge
                global nc

                self._atoms()
                charge = self.pahdb["species"][uid]["charge"]
                nc = self.atoms[uid]["nc"]

            if keywords.get('star') or keywords.get('isrf'):
                energy = Transitions.mean_energy(**keywords)
                self.model['energy'][uid]['sigma'] = np.sqrt(
                    Transitions.mean_energy_squared(**keywords) - energy**2
                )

            else:
                energy = e
                self.model['energy'][uid]['sigma'] = 0.

            self.model['energy'][uid]['e'] = energy

            print(
                f"MEAN ABSORBED ENERGY             : {self.model['energy'][uid]['e'] / 1.6021765e-12} +/- {self.model['energy'][uid]['sigma'] / 1.6021765e-12} eV"
            )

            global frequencies
            frequencies = np.array([d["frequency"] for d in self.data[uid]])

            if keywords.get('approximate'):
                Tmax = optimize.brentq(Transitions.approximate_attained_temperature, 2.73, 5000.0)
            else:
                Tmax = optimize.brentq(Transitions.attained_temperature, 2.73, 5000.0)

            self.model["temperatures"][uid]["temperature"] = Tmax

            print("MAXIMUM ATTAINED TEMPERATURE     : %f Kelvin" % Tmax)

            for d in self.data[uid]:
                if d["intensity"] > 0:
                    d["intensity"] *= (
                        2.4853427121856266e-23
                        * d["frequency"] ** 3
                        / (np.expm1(1.4387751297850830401 * d["frequency"] / Tmax))
                    )

            # Stop timer and calculate elapsed time.
            elapsed = timedelta(seconds=(time.perf_counter() - tstart))
            print(f"Elapsed time: {elapsed}")

            i += 1

        description = (
            f'model: calculated_temperature, approximated: {keywords.get("approximate", False)}'
        )

        if self.model['isrf']:
            description += ', isrf: yes'

        if self.model['star']:
            description += ', star: yes'
            description += f', Tstar: {Tstar} Kelvin'
            description += f', modelled: {keywords.get("stellar_model", False)}'

        else:
            description += f'<E>: {energy / 1.6021765e-12} eV'

        self.model['description'] = description

        print(57 * "=")

    def cascade(self, e: float, **keywords) -> None:
        """
        Applies the Cascade emission model.

        :param e: Excitation energy in erg.
        :type: float

        """
        if self.database != "theoretical" and not keywords.get('approximate'):
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
        TStar = 0.0

        tstart = time.perf_counter()

        if keywords.get('star') and keywords.get('stellarmodel'):
            message('STELLAR MODEL SELECTED: USING FIRST PARAMETER AS MODEL')

            TStar = (
                4
                * np.pi
                * integrate.simpson(energy["frequency"], energy["intensity"])
                / 5.67040e-5
            ) ** (0.25)

            select = np.where(
                (energy["frequency"] >= 2.5e3) & (energy["frequency"] <= 1.1e5)
            )
            nselect = energy["frequency"][select]

            if len(nselect) == 0:
                message('STELLAR MODEL HAS NO DATA BETWEEN 2.5E3-1.1E5 /cm')
                self.state = 0

            message('REBINNING STELLAR MODEL: 100 POINTS')
            # ^ Consider giving a user warning when providing large dimension data.

            global star_model
            star_model = [{'frequency': 0.0, 'intensity': 0}] * 100

            # * Using interpolation through ndimage.zoom to get equivalent
            # * of the IDL congrid function.
            star_model['frequency'] = ndimage.zoom(
                energy['frequency'][select], 100 / len(energy['frequency'][select])
            )
            star_model['intensity'] = ndimage.zoom(
                energy['intensity'][select], 100 / len(energy['intensity'][select])
            )

            message(f'CALCULATED EFFECTIVE TEMPERATURE: {TStar} Kelvin')

        elif keywords.get('star'):
            global Tstar
            Tstar = energy
            message(f'BLACKBODY TEMPERATURE: {TStar} Kelvin')

        if keywords.get('isrf') and not keywords.get('convolved'):
            message('ISRF SELECTED: IGNORING FIRST PARAMETER')

        if keywords.get('isrf') or keywords.get('star') and keywords.get('convolved'):
            message('CONVOLVING WITH ENTIRE RADIATION FIELD')

        self.model = {
            "type": "cascade_m",
            "energy": {},
            "temperatures": {},
            "approximate": keywords.get('approximate'),
            "star": keywords.get('star'),
            "isrf": keywords.get('isrf'),
            "convolved": keywords.get('convolved'),
            "stellar_model": keywords.get('stellar_model'),
            "Tstar": Tstar,
            "description": "",
        }

        self.units["ordinate"] = {"unit": u.erg, "label": "integrated radiant energy [erg]"}

        description = (
            f'model: cascade, approximated: {keywords.get("approximate", False)}'
        )

        if self.model['isrf']:
            description += ', isrf: yes'
            description += f', convolved: {keywords.get("convlolved", False)}'

        if self.model['star']:
            description += ', star: yes'
            description += f', Tstar: {Tstar} Kelvin'
            description += f', modelled: {keywords.get("stellar_model", False)}'
            description += f', convolved: {keywords.get("convolved", False)}'

        else:
            description += f'<E>: {energy / 1.6021765e-12} eV'

        self.model['description'] = description

        print(57 * "=")

        if keywords.get('approxiamte', False):
            message('USING APPROXIMATION')

            func1 = Transitions.approximate_attained_temperature
            func2 = Transitions.approximate_feature_strength

            if keywords.get('convolved', False):
                if keywords.get('isrf', False):
                    func3 = Transitions.isrf_approximate_feature_strength_convolved

                elif keywords.get('stellarmodel', False):
                    func3 = Transitions.stellar_model_approximate_feature_strength_convolved

                else:
                    func3 = Transitions.planck_approximate_feature_strength_convolved
        else:
            func1 = Transitions.attained_temperature
            func2 = Transitions.feature_strength

            if keywords.get('convolved', False):
                if keywords.get('isrf', False):
                    func3 = Transitions.isrf_feature_strength_convolved

                elif keywords.get('stellarmodel', False):
                    func3 = Transitions.stellar_model_feature_strength_convolved

                else:
                    func3 = Transitions.planck_feature_strength_convolved

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
                # Instatiate model energy and temperature dictionaries.
                self.model['energy'][uid] = {}
                self.model["temperatures"][uid] = {}

                print("SPECIES                          : %d/%d" % (i + 1, nuids))
                print("UID                              : %d" % uid)

                if keywords.get('approximate') or keywords.get('star') or keywords.get('isrf'):
                    global charge
                    global nc

                    self._atoms()
                    charge = [self.pahdb["species"][uid]["charge"] for uid in self.uids]
                    nc = np.array([self.atoms[uid]["nc"] for uid in self.uids])

                if keywords.get('star') or keywords.get('isrf'):
                    energy = Transitions.mean_energy(**keywords)

                    self.model['energy'][uid]['sigma'] = np.sqrt(
                        Transitions.mean_energy_squared(**keywords) - energy**2
                    )

                    if keywords.get('convolved'):
                        Nphot = Transitions.number_of_photons(**keywords)

                else:
                    energy = e
                    self.model['energy'][uid]['sigma'] = 0.

                self.model['energy'][uid]['e'] = energy

                print(
                    f"MEAN ABSORBED ENERGY             : {self.model['energy'][uid]['e'] / 1.6021765e-12} +/- {self.model['energy'][uid]['sigma'] / 1.6021765e-12} eV"
                )

                global frequencies
                frequencies = np.array([d["frequency"] for d in self.data[uid]])

                global intensities
                intensities = np.array([d["intensity"] for d in self.data[uid]])

                if keywords.get('approximate'):
                    totalcrossection = np.sum(intensities)

                Tmax = optimize.brentq(func1, 2.73, 5000.0)
                self.model["temperatures"][uid] = Tmax
                print("MAXIMUM ATTAINED TEMPERATURE     : %f Kelvin" % Tmax)

                if keywords.get('star') or keywords.get('isrf') and keywords.get('convolved'):
                    for d in self.data[uid]:
                        if d["intensity"] > 0:
                            global frequency
                            frequency = d["frequency"]
                            d["intensity"] *= (
                                d["frequency"] ** 3
                                * integrate.quad(func3, 2.5e3, 1.1e5)[
                                    0
                                ]
                            )
                    intensities = np.array([x / Nphot for x in intensities])

                else:
                    for d in self.data[uid]:
                        if d["intensity"] > 0:
                            global frequency
                            frequency = d["frequency"]
                            d["intensity"] *= (
                                d["frequency"] ** 3
                                * integrate.quad(func2, 2.73, Tmax)[
                                    0
                                ]
                            )

                i += 1

            if not keywords.get('approximate'):
                intensities = np.array([x * y**3 for x, y in zip(intensities, frequencies)])
            else:
                intensities = np.array([x * 2.48534271218563e-23 * nc * y**3 / totalcrossection for x, y in zip(intensities, frequencies)])

            print(f"ENERGY CONSERVATION IN SPECTRUM  :{np.sum(intensities) / self.model['energy'][uid]['e']}")

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
    def absorption_cross_section(f: np.ndarray) -> np.ndarray:
        """
        Calculates the PAH absorption cross-section per Li & Draine 2007,
        ApJ, 657:810-837.

        :Params f: frequencies in wavenumber

        :Returns: float array
        """

        wave = 1e4 / f

        A_c = [7.97e-17, 1.23e-17, 20e-21, 14e-21, 80e-24, 84e-24, 46e-24, -322e-24]
        W_c = [0.195, 0.217, 0.0805, 0.20, 0.0370, 0.0450, 0.0150, 0.135]
        C_c = [0.0722, 0.2175, 1.05, 1.23, 1.66, 1.745, 1.885, 1.90]

        A = np.transpose(np.resize(A_c, (np.size(f), np.size(A_c))))
        W = np.transpose(np.resize(W_c, (np.size(f), np.size(W_c))))
        C = np.transpose(np.resize(C_c, (np.size(f), np.size(C_c))))

        # Cutoff wavelength from Salama et al. (1996), over wavelength
        # (Eq. (4) in Mattioda et al. (2005))
        y = 1.0 / (0.889 + (2.282 / np.sqrt(0.4 * nc))) / wave

        wave_r2 = np.resize(wave, (2, (np.size(f))))

        crosssection = ((1.0 / np.pi) * np.arctan((1e3 * (y - 1.0) ** 3) / y) + 0.5) * (
            3458e-20 * 10.0 ** (-3.431 * wave)
            + (2.0 / np.pi)
            * np.sum(
                W[:2]
                * C[:2]
                * A[:2]
                / (((wave_r2 / C[:2]) - (C[:2] / wave_r2)) ** 2 + W[:2] ** 2),
                axis=0,
            )
        )

        if charge != 0:
            wave_r6 = np.resize(wave, (6, np.size(f)))

            crosssection = (
                crosssection
                + np.exp(-1e-1 / wave**2) * 1.5e-19 * 10 ** (-wave)
                + np.sqrt(2.0 / np.pi)
                * np.sum(
                    A[2:] * np.exp(-2.0 * (wave_r6 - C[2:]) ** 2 / W[2:] ** 2) / W[2:],
                    axis=0,
                )
            )

        return crosssection

    @staticmethod
    def planck(f: np.ndarray) -> np.ndarray:
        """
        Callback function to Calculate the absorption cross-section
        multiplied by Planck's function.

        :Params f: frequencies in wavenumber

        :Returns: float array

        """
        return (
            Transitions.absorption_cross_section(f)
            * f**3
            / (np.expm1(1.4387751297850830401 * f / Tstar))
        )

    @staticmethod
    def planck_squared(f: np.ndarray) -> np.ndarray:
        """
        Callback function to Calculate the absorption cross-section
        multiplied by Planck's function squared.

        :Params f: frequencies in wavenumber

        :Returns: float array

        """
        return (
            Transitions.absorption_cross_section(f)
            * f**4
            / (np.expm1(1.4387751297850830401 * f / Tstar))
        )

    @staticmethod
    def planck_number_of_photons(f: np.ndarray) -> np.ndarray:
        """
        Callback function to Calculate the number of photons
        using Planck's function.

        :Params f: frequencies in wavenumber

        :Returns: float array

        """
        return (
            Transitions.absorption_cross_section(f)
            * f**2
            / (np.expm1(1.4387751297850830401 * f / Tstar))
        )

    @staticmethod
    def isrf(f: np.ndarray) -> np.ndarray:
        """
        Callback function to calculate the absorption cross-section times
        the interstellar radiation field per Mathis et al. 1983, A&A, 128:212.

        :Params f: frequencies in wavenumber

        :Returns: float array
        """
        if f > 1.1e5:
            return 0.0
        if f > 1e4 / 0.110:
            return Transitions.absorption_cross_section(f) * 1.20e23 / f**5.4172
        if f > 1e4 / 0.134:
            return Transitions.absorption_cross_section(f) * 1.366e6 / f**2.0
        if f > 1e4 / 0.246:
            return Transitions.absorption_cross_section(f) * 1.019e-2 / f**0.3322

        T = [7500, 4000, 3000, 2.73]
        W = [1e-14, 1.65e-13, 4e-13, 1]

        return (
            Transitions.absorption_cross_section(f)
            * f**3
            * np.sum(
                [w / np.expm1(1.4387751297850830401 * f / t) for w, t in zip(W, T)]
            )
        )

    @staticmethod
    def isrf_squared(f: np.ndarray) -> np.ndarray:
        """
        Callback function to calculate the absorption cross-section times
        the interstellar radiation field per Mathis et al. 1983, A&A, 128:212.

        :Params f: frequencies in wavenumber

        :Returns: float array
        """
        if f > 1.1e5:
            return 0.0
        if f > 1e4 / 0.110:
            return Transitions.absorption_cross_section(f) * 1.20e23 / f**4.4172
        if f > 1e4 / 0.134:
            return Transitions.absorption_cross_section(f) * 1.366e6 / f
        if f > 1e4 / 0.246:
            return Transitions.absorption_cross_section(f) * 1.019e-2 / f**0.6678

        T = [7500, 4000, 3000, 2.73]
        W = [1e-14, 1.65e-13, 4e-13, 1]

        return (
            Transitions.absorption_cross_section(f)
            * f**4
            * np.sum(
                [w / np.expm1(1.4387751297850830401 * f / t) for w, t in zip(W, T)]
            )
        )

    @staticmethod
    def isrf_number_of_photons(f: np.ndarray) -> np.ndarray:
        """
        Callback function to calculate the number of photons per Mathis et al. 1983, A&A, 128:212.

        :Params f: frequencies in wavenumber

        :Returns: float array

        """
        if f > 1.1e5:
            return 0.0
        if f > 1e4 / 0.110:
            return Transitions.absorption_cross_section(f) * 1.20e23 / f**6.4172
        if f > 1e4 / 0.134:
            return Transitions.absorption_cross_section(f) * 1.366e6 / f**3
        if f > 1e4 / 0.246:
            return Transitions.absorption_cross_section(f) * 1.019e-2 / f**1.3322

        T = [7500, 4000, 3000, 2.73]
        W = [1e-14, 1.65e-13, 4e-13, 1]

        return (
            Transitions.absorption_cross_section(f)
            * f**2
            * np.sum(
                [w / np.expm1(1.4387751297850830401 * f / t) for w, t in zip(W, T)]
            )
        )

    @staticmethod
    def mean_energy(**keywords):
        """
        Callback function to calculate the mean energy in erg for a given
        blackbody temperatyre.

        :Returns: float array

        """
        # ! Using simpson's integration, given there is no Python (scipy)
        # ! implementation of the 5-point Newton-Cotes integration of
        # ! IDL's INT_TABULATED.

        if keywords.get('stellar_model'):
            me = (
                1.9864456023253396e-16
                * integrate.simpson(
                    star_model['frequency'],
                    Transitions.absorption_cross_section(star_model['frequency'])
                    * star_model['intensity'],
                )
                / integrate.simpson(
                    star_model['frequency'],
                    Transitions.absorption_cross_section(star_model['frequency'])
                    * star_model['intensity']
                    / star_model['frequency'],
                )
            )
        elif keywords.get('isrf'):
            me = (
                1.9864456023253396e-16
                * integrate.quad(Transitions.isrf, 2.5e3, 1.1e5)[0]
                / integrate.quad(Transitions.isrf_number_of_photons, 2.5e3, 1.1e5)[0]
            )
        else:
            me = (
                1.9864456023253396e-16
                * integrate.quad(Transitions.planck, 2.5e3, 1.1e5)[0]
                / integrate.quad(Transitions.planck_number_of_photons, 2.5e3, 1.1e5)[0]
            )

        return me

    @staticmethod
    def mean_energy_squared(**keywords):
        """
        Callback function to calculate the mean energy in erg^2 for a given
        blackbody temperatyre.

        :Returns: float array

        """
        # ! Using simpson's integration, given there is no Python (scipy)
        # ! implementation of the 5-point Newton-Cotes integration of
        # ! IDL's INT_TABULATED.

        if keywords.get('stellarmodel'):
            me = (
                3.945966130997681e-32
                * integrate.simpson(
                    star_model['frequency'],
                    star_model['frequency']
                    * Transitions.absorption_cross_section(star_model['frequency'])
                    * star_model['intensity'],
                )
                / integrate.simpson(
                    star_model['frequency'],
                    Transitions.absorption_cross_section(star_model['frequency'])
                    * star_model['intensity']
                    / star_model['frequency'],
                )
            )
        elif keywords.get('isrf'):
            me = (
                3.945966130997681e-32
                * integrate.quad(Transitions.isrf_squared, 2.5e3, 1.1e5)[0]
                / integrate.quad(Transitions.isrf_number_of_photons, 2.5e3, 1.1e5)[0]
            )
        else:
            me = (
                3.945966130997681e-32
                * integrate.quad(Transitions.planck_squared, 2.5e3, 1.1e5)[0]
                / integrate.quad(Transitions.planck_number_of_photons, 2.5e3, 1.1e5)[0]
            )

        return me

    @staticmethod
    def attained_temperature(T: float) -> float:
        """
        Calculate a PAH's temperature after absorbing a given amount of energy.

        :param T: Excitation temperature in Kelvin.
        :type T: float

        """
        global energy

        return integrate.quad(Transitions.heat_capacity, 2.73, T)[0] - energy

    @staticmethod
    def approximate_attained_temperature(T: float) -> float:
        """
        Calculate a PAH's temperature after absorbing a given amount of energy
        using an approximation.

        :param T: Excitation temperature in Kelvin.
        :type T: float

        """
        global energy

        return (
            nc
            * (
                7.54267e-11 * special.erf(-4.989231 + 0.41778 * np.log(T))
                + 7.542670e-11
            )
            - energy
        )

    @staticmethod
    def heat_capacity(T: float) -> float:
        """
        Calculate heat capacity.

        :param T: Excitation temperature in Kelvin.
        :type T: float

        """

        global frequencies

        val = 1.4387751297850830401 * frequencies / T

        return 1.3806505e-16 * np.sum(np.exp(-val) * (val / (1.0 - np.exp(-val))) ** 2)

    @staticmethod
    def planck_feature_strength_convolved(f: np.ndarray) -> np.ndarray:
        """
        Calculate a feature's strength convolved with a
        given blackbody radiation field.

        :Param f: frequencies in wavenumber

        :Returns: float

        """
        E = 1.9864456e-16 * f
        Tmax = optimize.brentq(Transitions.attained_temperature, 2.73, 5000.0)

        return (
            Transitions.planck_number_of_photons(f)
            * integrate.quad(Transitions.feature_strength, 2.73, Tmax)[0]
        )

    @staticmethod
    def isrf_feature_strength_convolved(f: np.ndarray) -> np.ndarray:
        """
        Calculate a feature's strength convolved with the
        interstellar radiation field.

        :Param f: frequencies in wavenumber

        :Returns: float

        """
        E = 1.9864456e-16 * f
        Tmax = optimize.brentq(Transitions.attained_temperature, 2.73, 5000.0)

        return (
            Transitions.isrf_number_of_photons(f)
            * integrate.quad(Transitions.feature_strength, 2.73, Tmax)[0]
        )

    @staticmethod
    def stellar_model_feature_strength_convolved(f: np.ndarray) -> np.ndarray:
        """
        Calculate a feature's strength convolved with a given stellar model.

        :Param f: frequencies in wavenumber

        :Returns: float

        """
        E = 1.9864456e-16 * f
        Tmax = optimize.brentq(Transitions.attained_temperature, 2.73, 5000.0)

        return (
            Transitions.absorption_cross_section(f)
            * np.interp(
                f,
                star_model['frequency'],
                star_model['intensity'] / star_model['frequenxy'],
            )
            * integrate.quad(Transitions.feature_strength, 2.73, Tmax)[0]
        )

    @staticmethod
    def planck_approximate_feature_strength_convolved(f: np.ndarray) -> np.ndarray:
        """
        Calculate a feature's strength convolved with a
        given blackbody using an approximation.

        :Param f: frequencies in wavenumber

        :Returns: float

        """
        E = 1.9864456e-16 * f
        Tmax = optimize.root_scalar(
            Transitions.attained_temperature, bracket=[2.73, 5000.0]
        )

        return (
            Transitions.planck_number_of_photons(f)
            * integrate.quad(Transitions.approximate_feature_strength, 2.73, Tmax)[0]
        )

    @staticmethod
    def isrf_approximate_feature_strength_convolved(f: np.ndarray) -> np.ndarray:
        """
        Calculate a feature's strength convolved with
        the interstellar radiation field using an approximation.

        :Param f: frequencies in wavenumber

        :Returns: float

        """
        E = 1.9864456e-16 * f
        Tmax = optimize.brentq(Transitions.attained_temperature, 2.73, 5000.0)

        return (
            Transitions.isrf_number_of_photons(f)
            * integrate.quad(Transitions.approximate_feature_strength, 2.73, Tmax)[0]
        )

    @staticmethod
    def stellar_model_approximate_feature_strength_convolved(
        f: np.ndarray,
    ) -> np.ndarray:
        """
        Calculate a feature's strength convolved with
        the interstellar radiation field using an approximation.

        :Param f: frequencies in wavenumber

        :Returns: float

        """
        E = 1.9864456e-16 * f
        Tmax = optimize.brentq(Transitions.attained_temperature, 2.73, 5000.0)

        return (
            Transitions.absorption_cross_section(f)
            * np.interp(
                f,
                star_model['frequency'],
                star_model['intensity'] / star_model['frequenxy'],
            )
            * integrate.quad(Transitions.approximate_feature_strength, 2.73, Tmax)[0]
            / 1.9864456023253396e-16
        )

    @staticmethod
    def feature_strength(T: float) -> float:
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

        return (Transitions.heat_capacity(T) / np.expm1(val1)) * (
            1.0
            / np.sum(
                intensities[valid] * (frequencies[valid]) ** 3 / np.expm1(val2[valid])
            )
        )

    @staticmethod
    def approximate_feature_strength(T: float):
        """
        Calculate a feature's strength convolved with a
        blackbody using an approximation from Bakes, Tielens & Bauschliher,
        ApJ, 556:501-514, 2001.

        : Param T: Excitation temperature in Kelvin.

        """
        a = 0.0
        b = 0.0

        if charge != 0:
            if T > 1000:
                a = 4.8e-4
                b = 1.6119

            elif (T > 300) and (T <= 1000):
                a = 6.38e-7
                b = 2.5556

            elif (T > 100) and (T <= 300):
                a = 1.69e-12
                b = 4.7687

            elif (T > 40) and (T <= 100):
                a = 7.7e-9
                b = 2.9244

            elif (T > 20) and (T <= 40):
                a = 3.4e-12
                b = 5.0428

            elif (T > 2.7) and (T <= 20):
                a = 4.47e-19
                b = 10.3870

        else:
            if T > 270:
                a = 5.5e-7
                b = 2.5270

            elif (T > 200) and (T <= 270):
                a = 1.7e-9
                b = 3.5607

            elif (T > 60) and (T <= 200):
                a = 1.35e-9
                b = 4.4800

            elif (T > 30) and (T <= 30):
                a = 4.18e-8
                b = 2.5217

            elif (T > 2.7) and (T <= 30):
                a = 1.8e-16
                b = 8.1860

        val = 1.4387751297850830401 * frequency / T

        if (a == 0.0) or (val > np.log(np.finfo(float).max)):
            return 0.0

        else:
            return 1 / ((np.exp(val) - 1) * a * T**b)

    @staticmethod
    def number_of_photons(**keywords):
        """
        Calculate the number of photons given a
        blacbody radiation field.

        """

        if keywords.get['stellar_model']:
            return integrate.simpson(
                star_model["frequency"],
                Transitions.absorption_cross_section(
                    star_model["frequency"]
                    * star_model["intensity"]
                    / star_model["frequency"]
                ),
            )

        elif keywords.get('isrf'):
            return integrate.quad(Transitions.isrf_number_of_photons, 2.5e3, 1.1e5)[0]

        else:
            return integrate.quad(Transitions.planck_number_of_photons, 2.5e3, 1.1e5)[0]

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
