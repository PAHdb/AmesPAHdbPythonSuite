#!/usr/bin/env python3
"""cluster_db_spectra.py

This is an example of clustering the spectra in the database built
around the functionality provided by the Suite and should help
confirm that the AmesPAHdbPythonSuite has been properly installed.

For more examples visit the PAHdb cookbook website:
https://pahdb.github.io/cookbook/

"""

import matplotlib.pyplot as plt
import numpy as np
from pkg_resources import resource_filename
from scipy.cluster.vq import kmeans
from scipy.integrate import simpson

from amespahdbpythonsuite.amespahdb import AmesPAHdb

if __name__ == "__main__":
    # Define FWHM, frequency range, number of points and number of clusters to use

    fwhm = 20.0

    npoints = 301

    xrange = 1e4 / np.array([15, 2.5])

    nclusters = 2

    # Read in the database.
    xml = "resources/pahdb-theoretical_cutdown.xml"
    pahdb = AmesPAHdb(
        filename=resource_filename("amespahdbpythonsuite", xml),
        check=False,
        cache=False,
    )

    # Retrieve all transitions for all species from the database
    uids = pahdb.search("c>20")

    transitions = pahdb.gettransitionsbyuid(uids)

    # Convolve spectra
    spectrum = transitions.convolve(fwhm=fwhm, xrange=xrange, npoints=npoints)

    # Obtain the spectrum for use outside object
    spectra = spectrum.get()

    #  Set all integrated intensities to unity
    for d in spectra["data"].values():
        d /= simpson(d, spectra["grid"])

    # Construct matrix
    matrix = np.array(list(spectra["data"].values()))

    # Perform k-means cluster analysis
    clusters, _ = kmeans(matrix, nclusters, iter=128)

    # Plot means
    for d in clusters:
        plt.plot(spectra["grid"], d)

    plt.xlabel("frequency [cm$^{\\rm -1}$]")
    plt.ylabel("normalized intensity [km/mol]")
    plt.show()
