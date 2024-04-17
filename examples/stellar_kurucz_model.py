#!/usr/bin/env python3
"""fit_a_spectrum.py

Program demonstrating fitting an astronomical spectrum using the
AmesPAHdbPythonSuite.

For more examples visit the PAHdb cookbook website:
https://pahdb.github.io/cookbook/

"""

import importlib_resources
import numpy as np

from amespahdbpythonsuite import observation
from amespahdbpythonsuite.amespahdb import AmesPAHdb

from astropy.io import fits


if __name__ == "__main__":
    file_path = importlib_resources.files("amespahdbpythonsuite")
    data_file = file_path / "resources/galaxy_spec.ipac"
    obs = observation.Observation(data_file)

    # Convert spectral unit to wavenumber required by AmesPAHdbPythonSuite.
    obs.abscissaunitsto("1/cm")

    # Read the database.
    xml_file = file_path / "resources/pahdb-theoretical_cutdown.xml"
    pahdb = AmesPAHdb(
        filename=xml_file,
        check=False,
        cache=False,
    )

    # Retrieve the transitions from the database for the subset of PAHs.
    uids = pahdb.search("c>0")
    transitions = pahdb.gettransitionsbyuid(uids)

    # Load the Kurucz stellar model
    fits_file = file_path / "resources/ckp00_17000.fits"
    with fits.open(fits_file) as hdulist:
        angstrom = hdulist[1].data["WAVELENGTH"]
        kurucz = {
            "frequency": 1e8 / np.flip(angstrom),
            "intensity": 1e-8
            * np.flip(hdulist[1].data["g40"] * angstrom**2)
            / (4 * np.pi),
        }

    # Calculate the emission spectrum convolved with a Kurucz stellar model.
    transitions.cascade(
        kurucz,
        star=True,
        stellar_model=True,
        convolved=True,
        multiprocessing=False,
        cache=False,
    )

    # Shift data 15 wavenumber to the red to mimic some effects of anharmonicity.
    transitions.shift(-15.0)

    # Convolve the bands with a Gaussian with FWHM of 15 /cm.
    spectrum = transitions.convolve(
        grid=obs.getgrid(), fwhm=15.0, gaussian=True, multiprocessing=False
    )

    # Fit the spectrum
    fit = spectrum.fit(obs)

    # Create plots.
    fit.plot(wavelength=True, residual=True)
    fit.plot(wavelength=True, size=True)
    fit.plot(wavelength=True, charge=True)
    fit.plot(wavelength=True, composition=True)

    # Predict 3 - 20 µm spectrum
    transitions.intersect(fit.getuids())

    xrange = 1e4 / np.array([20.0, 3.0])

    spectrum = transitions.convolve(xrange=xrange, gaussian=True, multiprocessing=False)

    coadded = spectrum.coadd(weights=fit.getweights())

    coadded.plot()
