#!/usr/bin/env python3
"""example.py

Example of using the AmesPAHdbPythonSuite to display the 'stick'
absorption spectrum of coronene (UID==18).

For more examples visit the PAHdb cookbook website:
https://pahdb.github.io/cookbook/

"""

import importlib_resources

from amespahdbpythonsuite.amespahdb import AmesPAHdb

if __name__ == "__main__":
    # Read the database.
    file_path = importlib_resources.files("amespahdbpythonsuite")
    xml_file = file_path / "resources/pahdb-theoretical_cutdown.xml"
    pahdb = AmesPAHdb(
        filename=xml_file,
        check=False,
        cache=False,
    )

    # Retrieve the transitions from the database for coronene.
    transitions = pahdb.gettransitionsbyuid([18])

    # Plot the emission 'stick' spectrum.
    transitions.plot(show=True)

    # Calculate the emission spectrum at the temperature reached
    # after absorbing a 6 eV (CGS units) photon.
    transitions.cascade(6 * 1.603e-12, multiprocessing=False)

    # Plot the emission 'stick' spectrum at that temperature.
    transitions.plot(show=True)

    # Convolve the bands with a Gaussian with FWHM of 15 /cm.
    convolved = transitions.convolve(fwhm=15.0, gaussian=True, multiprocessing=False)
    convolved.plot(show=True)
