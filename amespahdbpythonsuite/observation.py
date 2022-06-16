#!/usr/bin/env python3

import warnings
from typing import Union
import numpy as np

from astropy.io import fits, ascii

from astropy.io.registry import IORegistryError
from astropy.io.fits.verify import VerifyWarning
from astropy.nddata import StdDevUncertainty
from astropy import units as u
from specutils import Spectrum1D, SpectralRegion, manipulation

from amespahdbpythonsuite.amespahdb import AmesPAHdb

message = AmesPAHdb.message


class Observation:
    """
    AmesPAHdbPythonSuite observation class.
    Contains methods to work with astronomical spectra.

    """

    def __init__(self, d=None, **keywords):
        self.filepath = ""
        self.spectrum = Spectrum1D
        self.set(d, **keywords)

    def set(self, d, **keywords) -> None:
        """
        Populate properties.

        """

        if isinstance(d, str):
            self.read(d)
            return

        if d:
            if d.get("type", "") == self.__class__.__name__:
                if "filepath" not in keywords:
                    self.filepath = d["filepath"]
                if "spectrum" not in keywords:
                    self.spectrum = d["spectrum"]

        if "filepath" in keywords:
            self.filepath = keywords.get("filepath")
        if "spectrum" in keywords:
            self.spectrum = keywords.get("spectrum")

    def get(self) -> dict:
        """
        Assigns class variables to dictionary.

        """
        d = dict()
        d["type"] = self.__class__.__name__
        d["filepath"] = self.filepath
        d["spectrum"] = self.spectrum

        return d

    def read(self, filename: str) -> None:
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

            # Always work as if spectrum is a cube.
            if len(self.spectrum.flux.shape) == 1:
                self.spectrum = Spectrum1D(
                    flux=np.reshape(
                        self.spectrum.flux,
                        (
                            1,
                            1,
                        )
                        + self.spectrum.flux.shape,
                    ),
                    spectral_axis=self.spectrum.spectral_axis,
                )

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

                        # Create WCS object.
                        # self.wcs = wcs.WCS(hdu[0].header, naxis=2)

                        h0 = self.header["PS3_0"]
                        h1 = self.header["PS3_1"]

                        # Create Spectrum1D object.
                        flux = h.data.T * u.Unit(h.header["BUNIT"])
                        wave = hdu[h0].data[h1] * u.Unit(hdu[h0].columns[h1].unit)
                        self.spectrum = Spectrum1D(flux, spectral_axis=wave)

                        return None

                    # Use the WCS definitions for coordinate three
                    # linear.
                    if "CDELT3" in hdu_keys:

                        self.header = h.header

                        # Create WCS object
                        # self.wcs = wcs.WCS(hdu[0].header, naxis=2)

                        # Create Spectrum1D object
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
            # Create Spectrum1D object.
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
        raise OSError(self.filepath + ": Format not recognized")

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
            g = x

        resampler = manipulation.FluxConservingResampler(
            extrapolation_treatment="nan_fill"
        )

        self.spectrum = resampler(self.spectrum, g * self.spectrum.spectral_axis.unit)

    def setgridrange(self, min: float, max=None) -> None:
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
