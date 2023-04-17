#!/usr/bin/env python3

from __future__ import annotations

from typing import TYPE_CHECKING, Optional, Union

import astropy.units as u  # type: ignore
import numpy as np
from astropy.nddata import StdDevUncertainty  # type: ignore
from scipy import optimize  # type: ignore
from specutils import Spectrum1D, manipulation  # type: ignore

from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite.transitions import Transitions

if TYPE_CHECKING:
    from amespahdbpythonsuite.coadded import Coadded
    from amespahdbpythonsuite.fitted import Fitted
    from amespahdbpythonsuite.observation import Observation
    from amespahdbpythonsuite.mcfitted import MCFitted

import multiprocessing as mp
from functools import partial

message = AmesPAHdb.message


class Spectrum(Transitions):
    """
    AmesPAHdbPythonSuite spectrum class.
    Contains methods to fit and plot the input spectrum.

    """

    def __init__(self, d: Optional[dict] = None, **keywords) -> None:
        super().__init__(d, **keywords)
        self.__set(d, **keywords)

    def set(self, d: Optional[dict] = None, **keywords) -> None:
        """
        Calls class: :class:`amespahdbpythonsuite.transitions.Transitions.set` to parse keywords.

        """
        Transitions.set(self, d, **keywords)
        self.__set(d, **keywords)

    def __set(self, d: Optional[dict] = None, **keywords) -> None:
        """
        Populate data dictionary helper.

        """
        self.grid = keywords.get("grid", list())
        self.profile = keywords.get("profile", "")
        self.fwhm = keywords.get("fwhm", 0.0)

        if isinstance(d, dict):
            if d.get("type", "") == self.__class__.__name__:
                if "grid" not in keywords:
                    self.grid = d["grid"]
                if "profile" not in keywords:
                    self.profile = d["profile"]
                if "fwhm" not in keywords:
                    self.fwhm = d["fwhm"]

    def get(self) -> dict:
        """
        Calls class: :class:`amespahdbpythonsuite.transitions.Transitions.get`.
        Assigns class variables from inherited dictionary.

        """
        d = Transitions.get(self)
        d["type"] = self.__class__.__name__
        d["grid"] = self.grid
        d["profile"] = self.profile
        d["fwhm"] = self.fwhm

        return d

    def __repr__(self) -> str:
        """
        Class representation.

        """
        return (
            f"{self.__class__.__name__}("
            f"{self.uids=},{self.grid=},{self.profile=},{self.fwhm=})"
        )

    def __str__(self) -> str:
        """
        A description of the instance.
        """

        return f"AmesPAHdbPythonSuite Spectrum instance.\n" f"{self.uids=}"

    def write(self, filename: str = "") -> None:
        """
        Write the spectra to file as an IPAC-table.

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
                np.array([f for _ in self.data.values() for f in self.grid])
                * self.units["abscissa"]["unit"],
                np.array([t for v in self.data.values() for t in v])
                * self.units["ordinate"]["unit"],
            ],
            names=["UID", "FREQUENCY", "INTENSITY"],
            meta={"comments": hdr},
        )

        ascii.write(tbl, filename, format="ipac", overwrite=True)

        message(f"WRITTEN: {filename}")

    def fit(
        self,
        y: Union[Observation, Spectrum1D, list],
        yerr: list = list(),
        notice: bool = True,
        **keywords,
    ) -> Optional[Fitted]:
        """
        Fits the input spectrum.

        """

        from amespahdbpythonsuite import observation

        if isinstance(y, Spectrum1D):
            obs = y
        elif isinstance(y, observation.Observation):
            obs = y.spectrum
        else:
            unc = None
            if np.any(yerr):
                unc = StdDevUncertainty(yerr)
            obs = Spectrum1D(
                flux=y * u.Unit(),
                spectral_axis=self.grid * self.units["abscissa"]["unit"],
                uncertainty=unc,
            )

        if obs.spectral_axis.unit != u.Unit() and obs.spectral_axis.unit != u.Unit(
            "1/cm"
        ):
            message("EXPECTING SPECTRAL UNITS OF 1 / CM")
            return None

        matrix = np.array(list(self.data.values()))

        if obs.uncertainty is None:
            method = "NNLS"
            b = list(obs.flux.value)
            m = matrix
        else:
            method = "NNLC"
            b = list(np.divide(obs.flux.value, obs.uncertainty.array))
            m = np.divide(matrix, obs.uncertainty.array)

        if notice:
            message(f"DOING {method}")

        solution, _ = optimize.nnls(m.T, b)

        # Initialize lists and dictionaries.
        uids = list()
        data = dict()
        weights = dict()

        # Retrieve uids, data, and fit weights dictionaries.
        for (
            uid,
            s,
            m,
        ) in zip(self.uids, solution, matrix):
            if s > 0:
                intensities = []
                uids.append(uid)
                for d in m:
                    intensities.append(s * d)
                data[uid] = np.array(intensities)
                weights[uid] = s

        if notice:
            message(
                [
                    " NOTICE: PLEASE TAKE CONSIDERABLE CARE WHEN INTERPRETING ",
                    " THESE RESULTS AND PUTTING THEM IN AN ASTRONOMICAL       ",
                    " CONTEXT. THERE ARE MANY SUBTLETIES THAT NEED TO BE TAKEN",
                    " INTO ACCOUNT, RANGING FROM PAH SIZE, INCLUSION OF       ",
                    " HETEROATOMS, ETC. TO DETAILS OF THE APPLIED EMISSION    ",
                    " MODEL, BEFORE ANY THOROUGH ASSESSMENT CAN BE MADE.      ",
                ]
            )

        from amespahdbpythonsuite.fitted import Fitted

        return Fitted(
            database=self.database,
            version=self.version,
            data=data,
            pahdb=self.pahdb,
            uids=uids,
            model=self.model,
            units=self.units,
            shift=self._shift,
            grid=self.grid,
            profile=self.profile,
            fwhm=self.fwhm,
            observation=obs,
            weights=weights,
            method=method,
        )

    def plot(self, **keywords) -> None:
        """
        Plot the spectrum.

        """
        import matplotlib.cm as cm  # type: ignore
        import matplotlib.pyplot as plt  # type: ignore

        _, ax = plt.subplots()
        ax.minorticks_on()
        ax.tick_params(which="major", right="on", top="on", direction="in")
        colors = cm.rainbow(np.linspace(0, 1, len(self.uids)))
        for y, col in zip(self.data.values(), colors):
            ax.plot(self.grid, y, color=col)

        ax.set_xlim((max(self.grid), min(self.grid)))

        ax.set_xlabel(
            self.units["abscissa"]["label"]
            + " ["
            + self.units["abscissa"]["unit"].to_string("latex_inline")
            + "]",
        )

        unit = self.units["ordinate"]["unit"]
        scale = unit.scale
        unit /= scale
        unit = unit.decompose().cgs.unit
        pre = ""

        if scale != 1.0:
            s = np.log10(scale)
            pre = r"$\times10^{" + f"{s:.0f}" + r"}$ "

        ax.set_ylabel(
            self.units["ordinate"]["label"]
            + " ["
            + pre
            + unit.to_string("latex_inline")
            + "]",
        )

        basename = keywords.get("save")
        if basename:
            if not isinstance(basename, str):
                basename = "spectrum"
            plt.savefig(f"{basename}.pdf")
        elif keywords.get("show", False):
            plt.show()

    def coadd(self, weights: dict = dict(), average: bool = False) -> Coadded:
        """
        Co-add PAHdb spectra.

        Parameters:
            weights: dict
                Dictionary of fit weights to use when co-adding.
            average: bool
                If True calculates the average coadded spectrum.

        """

        data: Union[np.ndarray, float]

        if weights:
            data = np.zeros(len(self.grid))
            for uid, weight in weights.items():
                data += self.data[uid] * weight
        else:
            data = sum(self.data.values())

        if average:
            data /= len(self.data.keys())

        from amespahdbpythonsuite.coadded import Coadded

        return Coadded(
            database=self.database,
            version=self.version,
            data={0: data},
            pahdb=self.pahdb,
            uids=[0],
            model=self.model,
            units=self.units,
            shift=self._shift,
            grid=self.grid,
            profile=self.profile,
            fwhm=self.fwhm,
            weights=weights,
            averaged=average,
        )

    def normalize(self, all: bool = False) -> Union[float, dict]:
        """
        Normalize spectral data

        """

        max: Union[float, dict] = 0.0
        if all:
            for intensities in self.data.values():
                m = intensities.max()
                if m > max:
                    max = m
            for intensities in self.data.values():
                intensities /= max
        else:
            max = dict()
            for uid, intensities in self.data.items():
                m = intensities.max()
                intensities /= m
                max[uid] = m

        return max

    def resample(self, grid: np.ndarray) -> None:
        """
        Resample the spectral data.

        """

        resampler = manipulation.FluxConservingResampler(
            extrapolation_treatment="nan_fill"
        )
        for uid, intensities in self.data.items():
            s = resampler(
                Spectrum1D(
                    spectral_axis=self.grid * self.units["abscissa"]["unit"],
                    flux=intensities * self.units["ordinate"]["unit"],
                ),
                grid * self.units["abscissa"]["unit"],
            )
            self.data[uid] = s.flux.value

        self.grid = grid

    def mcfit(
        self,
        y: Union[Observation, Spectrum1D, list],
        yerr: list = list(),
        samples: int = 1024,
        uniform: bool = False,
        multiprocessing: bool = False,
        notice: bool = True,
        **keywords,
    ) -> Optional[MCFitted]:
        """
        Monte Carlo sampling and fitting to the input spectrum.

        Parameters
        ----------
        samples : Number of samples.
            int

        """
        from tqdm import tqdm  # type: ignore

        from amespahdbpythonsuite import observation
        from amespahdbpythonsuite.mcfitted import MCFitted

        if isinstance(y, Spectrum1D):
            obs = y
        elif isinstance(y, observation.Observation):
            obs = y.spectrum
        else:
            unc = None
            if np.any(yerr):
                unc = StdDevUncertainty(yerr)
            obs = Spectrum1D(
                flux=y * u.Unit(),
                spectral_axis=self.grid * self.units["abscissa"]["unit"],
                uncertainty=unc,
            )

        if obs.spectral_axis.unit != u.Unit() and obs.spectral_axis.unit != u.Unit(
            "1/cm"
        ):
            message("EXPECTING SPECTRAL UNITS OF 1 / CM")
            return None

        if obs.uncertainty is None:
            message("UNCERTAINTIES REQUIRED FOR MCFIT")
            return None

        if notice:
            message(
                [
                    " NOTICE: PLEASE TAKE CONSIDERABLE CARE WHEN INTERPRETING ",
                    " THESE RESULTS AND PUTTING THEM IN AN ASTRONOMICAL       ",
                    " CONTEXT. THERE ARE MANY SUBTLETIES THAT NEED TO BE TAKEN",
                    " INTO ACCOUNT, RANGING FROM PAH SIZE, INCLUSION OF       ",
                    " HETEROATOMS, ETC. TO DETAILS OF THE APPLIED EMISSION    ",
                    " MODEL, BEFORE ANY THOROUGH ASSESSMENT CAN BE MADE.      ",
                ]
            )

        mcfits = list()

        # Start the MC sampling and fitting.
        if multiprocessing:
            from amespahdbpythonsuite.fitted import Fitted

            matrix = np.array(list(self.data.values()))

            m = np.divide(matrix, obs.uncertainty.array)

            pool = mp.Pool(mp.cpu_count() - 1)
            for solution, b in tqdm(
                pool.imap_unordered(
                    partial(
                        _mcfit,
                        m=m,
                        x=obs.flux.value,
                        u=obs.uncertainty.array,
                        uniform=uniform,
                    ),
                    range(samples),
                ),
                desc="samples",
                leave=True,
                unit="samples",
                colour="blue",
                total=samples,
            ):
                # Initialize lists and dictionaries.
                uids = list()
                data = dict()
                weights = dict()

                # Retrieve uids, data, and fit weights dictionaries.
                for (
                    uid,
                    s,
                    m,
                ) in zip(self.uids, solution, matrix):
                    if s > 0:
                        intensities = []
                        uids.append(uid)
                        for d in m:
                            intensities.append(s * d)
                        data[uid] = np.array(intensities) * obs.flux.unit
                        weights[uid] = s

                obs_fit = Spectrum1D(
                    flux=b * obs.flux.unit,
                    spectral_axis=obs.spectral_axis,
                    uncertainty=obs.uncertainty,
                )

                mcfits.append(
                    Fitted(
                        database=self.database,
                        version=self.version,
                        data=data,
                        pahdb=self.pahdb,
                        uids=uids,
                        model=self.model,
                        units=self.units,
                        shift=self._shift,
                        grid=self.grid,
                        profile=self.profile,
                        fwhm=self.fwhm,
                        observation=obs_fit,
                        weights=weights,
                        method="NNLC",
                    )
                )
        else:
            for _ in tqdm(
                range(samples),
                desc="samples",
                leave=True,
                unit="samples",
                colour="blue",
            ):
                if uniform:
                    # Calculate new flux based on random uniform distribution sampling.
                    flux = (
                        obs.uncertainty.array * np.random.uniform(-1, 1, obs.flux.shape)
                        + obs.flux
                    )
                else:
                    # Calculate new flux based on random normal distribution sampling.
                    flux = np.random.normal(obs.flux.value, obs.uncertainty.array)

                # Fit the spectrum.
                fit = self.fit(flux * obs.flux.unit, obs.uncertainty, notice=False)

                # Obtain the fit and weights.
                if fit:
                    mcfits.append(fit)

        return MCFitted(
            mcfits=mcfits,
            distribution="uniform" if uniform else "normal",
            observation=obs,
        )


def _mcfit(_, m, x, u, uniform) -> tuple:
    if uniform:
        # Calculate new flux based on random uniform distribution sampling.
        b = u * np.random.uniform(-1, 1, x.shape) + x
    else:
        # Calculate new flux based on random normal distribution sampling.
        b = np.random.normal(x, u)
    b = list(np.divide(b, u))

    # Fit the spectrum.
    solution, _ = optimize.nnls(m.T, b)

    return solution, b
