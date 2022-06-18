#!/usr/bin/env python3

from typing import Optional, Union

import numpy as np

from specutils import Spectrum1D, manipulation  # type: ignore
from astropy.nddata import StdDevUncertainty  # type: ignore
import astropy.units as u  # type: ignore
from scipy import optimize  # type: ignore

from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite.transitions import Transitions

message = AmesPAHdb.message


class Spectrum(Transitions):
    """
    AmesPAHdbPythonSuite spectrum class.
    Contains methods to fit and plot the input spectrum.

    """

    def __init__(self, d: Optional[dict] = None, **keywords) -> None:
        super().__init__(d, **keywords)
        self.__set(d, **keywords)

    def set(self, d, **keywords):
        """
        Calls class: :class:`amespahdbpythonsuite.transitions.Transitions.set` to parse keywords.

        """
        Transitions.set(self, d, **keywords)
        self.__set(d, **keywords)

    def __set(self, d: Optional[dict] = None, **keywords) -> None:
        """
        Populate data dictionary helper.

        """
        if isinstance(d, dict):
            if d.get("type", "") == self.__class__.__name__:
                if "grid" not in keywords:
                    self.grid = d["grid"]
                if "profile" not in keywords:
                    self.profile = d["profile"]
                if "fwhm" not in keywords:
                    self.fwhm = d["fwhm"]

        self.grid = keywords.get("grid", list())
        self.profile = keywords.get("profile", "")
        self.fwhm = keywords.get("fwhm", 0.0)

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

    def fit(self, y: list, yerr: list = [], notice: bool = True, **keywords):
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

        obs.spectral_axis.to("1/cm")

        matrix = np.array(list(self.data.values()))

        if obs.uncertainty is None:
            method = "NNLS"
            b = list(obs.flux.value)
            m = matrix
        else:
            method = "NNLC"
            b = list(np.divide(obs.flux.value, obs.uncertainty.array))
            m = np.divide(matrix, obs.uncertainty.array)

        message(f"DOING {method}")

        solution, norm = optimize.nnls(m.T, b)

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
            type=self.type,
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
        import matplotlib.pyplot as plt  # type: ignore
        import matplotlib.cm as cm  # type: ignore

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

    def coadd(self, weights: list = [], average: bool = False):
        """
        Co-add PAHdb spectra.

        Parameters:
            weights: dict
                Dictionary of fit weights to use when co-adding.
            average: bool
                If True calculates the average coadded spectrum.

        """

        data = np.zeros(len(self.grid))

        if weights:
            for key in self.data.keys():
                data += self.data[key] * weights[key]
        else:
            for key in self.data.keys():
                data += self.data[key]

        if average:
            data /= len(self.data.keys())

        from amespahdbpythonsuite.coadded import Coadded

        return Coadded(
            type=self.type,
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

        max = 0.0  # type: Union[float, dict]
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
