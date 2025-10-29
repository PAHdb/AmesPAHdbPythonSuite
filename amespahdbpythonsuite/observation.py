#!/usr/bin/env python3

import warnings
from pathlib import Path
from typing import Optional, Union

import numpy as np
from astropy import units as u  # type: ignore
from astropy.io import ascii, fits  # type: ignore
from astropy.io.fits.verify import VerifyWarning  # type: ignore
from astropy.io.registry import IORegistryError  # type: ignore
from astropy.nddata import StdDevUncertainty  # type: ignore
from specutils import SpectralRegion, Spectrum1D, manipulation  # type: ignore

from amespahdbpythonsuite.amespahdb import AmesPAHdb

message = AmesPAHdb.message


class Observation:
    """
    AmesPAHdbPythonSuite observation class.
    Contains methods to work with astronomical spectra.

    """

    filepath: Union[str, Path] = ""
    spectrum = Spectrum1D

    def __init__(self, d: Optional[None] = None, **keywords) -> None:
        self.set(d, **keywords)

    def set(self, d, **keywords) -> None:
        """
        Populate properties.

        """

        if isinstance(d, (Path, str)):
            self.read(d)
            return

        if isinstance(d, dict):
            if d.get("type", "") == self.__class__.__name__:
                if "filepath" not in keywords:
                    self.filepath = d["filepath"]
                if "spectrum" not in keywords:
                    self.spectrum = d["spectrum"]

        filepath = keywords.get("filepath")
        if "filepath" and isinstance(filepath, (Path, str)):
            self.filepath = filepath
        if "spectrum" in keywords:
            self.spectrum = keywords.get("spectrum")

    def get(self) -> dict:
        """
        Assigns class variables to dictionary.

        """
        d: dict[str, Union[str, Path, Spectrum1D]] = dict()
        d["type"] = self.__class__.__name__
        d["filepath"] = self.filepath
        d["spectrum"] = self.spectrum

        return d

    def __repr__(self) -> str:
        """
        Class representation.

        """
        return f"{self.__class__.__name__}(" f"{self.filepath=})"

    def __str__(self) -> str:
        """
        A description of the instance.
        """

        return f"AmesPAHdbPythonSuite Observation instance.\n" f"{self.filepath=}"

    def write(self, filename: str = "") -> None:
        """
        Write the spectrum to file as an IPAC-table.

        """
        import datetime
        import sys

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
        }

        for key, value in kv.items():
            if not value.isnumeric():
                hdr.append(f"{key:8} = '{value}'")
            else:
                hdr.append(f"{key:8} = {value}")

        tbl = Table(
            [
                self.spectrum.spectral_axis.quantity,
                self.spectrum.flux,
            ],
            names=["ABSCISSA", "ORDINATE"],
            meta={"comments": hdr},
        )

        if self.spectrum.uncertainty:
            tbl.add_column(self.spectrum.uncertainty.quantity, name="UNCERTAINTY")

        ascii.write(tbl, filename, format="ipac", overwrite=True)

        message(f"WRITTEN: {filename}")

    def plot(self, **keywords) -> None:
        """
        Plot the spectrum.
        """

        import matplotlib.pyplot as plt  # type: ignore

        if self.spectrum.uncertainty:
            plt.errorbar(
                self.spectrum.spectral_axis.value,
                self.spectrum.flux.value,
                self.spectrum.uncertainty.array,
                capsize=2,
            )
        else:
            plt.plot(self.spectrum.spectral_axis, self.spectrum.flux)
        plt.xlabel(self.spectrum.spectral_axis.unit.to_string("latex_inline"))
        plt.ylabel(self.spectrum.flux.unit.to_string("latex_inline"))

        basename = keywords.get("save")
        if basename:
            if not isinstance(basename, str):
                basename = "laboratory"
            plt.savefig(f"{basename}.pdf")
        elif keywords.get("show", False):
            plt.show()

    def read(self, filename: Union[Path, str]) -> None:
        """
        Read a spectrum.

        Parameters:
          filename: str
              Name of file to read.

        """

        self.filepath = filename

        try:
            # Supress warning when Spectrum1D cannot load the file.
            warnings.simplefilter("ignore", category=VerifyWarning)

            self.spectrum = Spectrum1D.read(self.filepath)

            if "header" in self.spectrum.meta:
                self.header = self.spectrum.meta["header"]
            else:
                self.header = fits.header.Header()

            return None
        except FileNotFoundError as e:
            raise (e)
        except (OSError, IORegistryError):
            # Because Spectrum1D raises a generic OSError when the
            # file cannot be read, we have to catch OSError here and pass
            # so that we can try and read it directly as FITS or ASCII.
            pass

        try:
            with fits.open(self.filepath) as hdu:
                for h in hdu:
                    hdu_keys = list(h.header.keys())

                    # Use the WCS definitions for coordinate three
                    # lookup table.
                    if "PS3_0" in hdu_keys and "PS3_1" in hdu_keys:
                        self.header = h.header

                        # Create WCS instance.
                        # self.wcs = wcs.WCS(hdu[0].header, naxis=2)

                        h0 = self.header["PS3_0"]
                        h1 = self.header["PS3_1"]

                        # Create Spectrum1D instance.
                        flux = h.data.T * u.Unit(h.header["BUNIT"])
                        wave = hdu[h0].data[h1].squeeze() * u.Unit(
                            hdu[h0].columns[h1].unit
                        )
                        self.spectrum = Spectrum1D(flux, spectral_axis=wave)

                        return None

                    # Use the WCS definitions for coordinate three
                    # linear.
                    if "CDELT3" in hdu_keys:
                        self.header = h.header

                        # Create WCS instance
                        # self.wcs = wcs.WCS(hdu[0].header, naxis=2)

                        # Create Spectrum1D instance
                        # u.Unit(self.header['BUNIT'])
                        flux = h.data.T * u.Unit("Jy")
                        wave = (
                            h.header["CRVAL3"]
                            + h.header["CDELT3"] * np.arange(0, h.header["NAXIS3"])
                        ) * u.Unit(h.header["CUNIT3"])
                        self.spectrum = Spectrum1D(flux, spectral_axis=wave)

                        return None

        except OSError:
            # Because astropy.io.fits.open raises a generic OSError
            # when the file header is missing the END card (which
            # ASCII files do), we have to catch OSError here and pass
            # so that we can try and read it as ASCII.
            pass

        try:
            data = ascii.read(self.filepath)
            for name in data.colnames:
                data.rename_column(name, name.upper())
            unc = None
            if "FLUX_UNCERTAINTY" in data.colnames:
                unc = StdDevUncertainty(data["FLUX_UNCERTAINTY"].quantity)
            # Create Spectrum1D instance.
            self.spectrum = Spectrum1D(
                flux=data["FLUX"].quantity,
                spectral_axis=data["WAVELENGTH"].quantity,
                uncertainty=unc,
            )
            str = ""
            for card in data.meta["keywords"].keys():
                value = data.meta["keywords"][card]["value"]
                str += "%-8s=%71s" % (card, value)
            self.header = fits.header.Header.fromstring(str)
            return None
        except Exception:
            pass

        # Like astropy.io we, simply raise a generic OSError when
        # we fail to read the file.
        raise OSError(f"{self.filepath}: Format not recognized")

    def getgrid(self) -> np.ndarray:
        """
        Retreive the spectral axis.

        """

        return self.spectrum.spectral_axis.value

    def rebin(
        self, x: Union[np.ndarray, float], uniform=False, resolution=False
    ) -> None:
        """
        Resample the spectral data.

        """

        g: Union[np.ndarray, list, float] = x
        if uniform or resolution:
            min = self.spectrum.spectral_axis.value.min()
            max = self.spectrum.spectral_axis.value.max()
            if uniform:
                message(f"REBINNING TO UNIFORM GRID: DELTA={x}")
                g = np.arange(min, max, x)
                if g[-1] != max:
                    g = np.append(g, max)
            elif resolution:
                message(f"REBINNING TO RESOLUTION: R={x}")
                g = [min]
                while g[-1] < max:
                    g.append(g[-1] + g[-1] / x)
                g[-1] = max
                g = np.array(g)
        else:
            message("REBINNING TO SET GRID")

        resampler = manipulation.FluxConservingResampler(
            extrapolation_treatment="nan_fill"
        )

        self.spectrum = resampler(self.spectrum, g * self.spectrum.spectral_axis.unit)

    def setgridrange(self, min: float, max: Optional[float] = None) -> None:
        """
        Truncate the data to the given range.

        """

        if not max:
            max = self.spectrum._spectral_axis.value.max()

        u = self.spectrum._spectral_axis.unit

        self.spectrum = manipulation.extract_region(
            self.spectrum, SpectralRegion(min * u, max * u)
        )

    def abscissaunitsto(self, unit: u.Unit) -> None:
        """
        Convert abscissa units.

        """

        self.spectrum = Spectrum1D(
            flux=self.spectrum.flux,
            spectral_axis=self.spectrum.spectral_axis.to(unit),
            uncertainty=self.spectrum.uncertainty,
        )
