#!/usr/bin/env python3
"""fit_a_spectrum.py

Program demonstrating using a stellar Kurucz model as the radiation source when
computing a PAH emission spectrum with the AmesPAHdbPythonSuite.

For more examples visit the PAHdb cookbook website:
https://pahdb.github.io/cookbook/

"""

import importlib_resources
import numpy as np

from amespahdbpythonsuite.amespahdb import AmesPAHdb

from astropy.io import fits


if __name__ == "__main__":

    # Read the database.
    file_path = importlib_resources.files("amespahdbpythonsuite")
    xml_file = file_path / "resources/pahdb-theoretical_cutdown.xml"
    pahdb = AmesPAHdb(
        filename=xml_file,
        check=False,
        cache=False,
    )

    # Retrieve the transitions from the database for the subset of PAHs.
    uid = 18

    transitions = pahdb.gettransitionsbyuid(uid)

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
    )

    # Shift data 15 wavenumber to the red to mimic some effects of anharmonicity.
    transitions.shift(-15.0)

    # Convolve the bands with a Gaussian with FWHM of 15 /cm.
    spectrum = transitions.convolve(fwhm=15.0, lorentzian=True, multiprocessing=False)

    # Plot the spectrum
    spectrum.plot(show=True, legend=True)
