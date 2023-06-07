#!/usr/bin/env python3
"""cluster_db_spectra.py

This is an example of producing a number of random-mixture PAH
emission spectrum, built around the functionality provided by the
AmesPAHdbPythonSuite and should help confirm that the it has been
properly installed.

For more examples visit the PAHdb cookbook website:
https://pahdb.github.io/cookbook/

"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from amespahdbpythonsuite.amespahdb import AmesPAHdb

if __name__ == "__main__":
    # Test parameters
    n = 1000

    # Spectral parameters
    Ein = 6.0 * 1.6021765e-12

    fwhm = 15.0

    gaussian = 0

    wrange = np.array([2.5, 20])

    npoints = 300

    # retrieve spectra from the database
    xml = "resources/pahdb-theoretical_cutdown.xml"
    pahdb = AmesPAHdb(
        #   filename=resource_filename("amespahdbpythonsuite", xml),
        check=False,
        #    cache=False,
    )

    # Retrieve all transitions for all species from the database
    uids = pahdb.search("mg=0 fe=0 si=0 o=0 c>20 chx=0 ch2=0")

    transitions = pahdb.gettransitionsbyuid(uids)

    transitions.cascade(Ein, multiprocessing=False)

    spectrum = transitions.convolve(fwhm=fwhm, xrange=1e4 / wrange, npoints=npoints)

    wavelength = 1e4 / spectrum.getgrid()

    spectra = spectrum.get()

    species = pahdb.getspeciesbyuid(uids)

    properties = species.get()

    # run test
    samples = np.zeros((n, npoints))

    sizeC = np.zeros(n)

    sizeH = np.zeros(n)

    charge = np.zeros((5, n))

    for i in range(n):
        print(f"\r{i+1:4d}", end="")

        abundances = np.random.uniform(size=len(uids))

        for abundance, uid in zip(abundances, uids):
            samples[i, :] += abundance * spectra["data"][uid]

            sizeC[i] += abundance * properties["data"][uid]["n_c"]

            sizeH[i] += abundance * properties["data"][uid]["n_h"]

            charge[properties["data"][uid]["charge"] + 1, i] += abundance

    matrix = np.zeros((n, n))
    for i in range(0, n):
        for j in range(i, n):
            if i == j:
                matrix[i, j] = 1.0
                continue
            matrix[i, j] = stats.pearsonr(samples[i, :], samples[j, :]).statistic
            matrix[j, i] = matrix[i, j]

    print("\rdone!")

    for sample in samples:
        sample /= sum(sample)

    moments = stats.describe(samples, axis=0)

    h, locations = np.histogram(matrix, bins="auto")

    imax = np.argmax(h)

    select = np.where(matrix != 1)

    avg = stats.describe(matrix[select])

    med = np.median(matrix[select])

    hC, locationsC = np.histogram(sizeC / len(uids), bins="auto")

    imaxC = np.argmax(hC)

    avgC = stats.describe(sizeC / len(uids))

    hH, locationsH = np.histogram(sizeH / len(uids), bins="auto")

    imaxH = np.argmax(hH)

    avgH = stats.describe(sizeH / len(uids))

    hCharge, locationsCharge = np.histogram(
        np.sum(charge[2:, :], axis=0) / np.sum(charge, axis=0), bins="auto"
    )

    imaxCharge = np.argmax(hCharge)

    avgCharge = stats.describe(np.sum(charge[2:, :], axis=0) / np.sum(charge, axis=0))

    hDbC, locationsDbC = np.histogram(
        [d["n_c"] for d in properties["data"].values()], bins="auto"
    )

    imaxDbC = np.argmax(hDbC)

    avgDbC = stats.describe([d["n_c"] for d in properties["data"].values()])

    hDbH, locationsDbH = np.histogram(
        [d["n_h"] for d in properties["data"].values()], bins="auto"
    )

    imaxDbh = np.argmax(hDbC)

    avgDbH = stats.describe([d["n_h"] for d in properties["data"].values()])

    # Plot results
    for sample in samples:
        plt.plot(wavelength, sample, color="grey")

    plt.fill_between(
        wavelength,
        moments.mean - np.sqrt(moments.variance),
        moments.mean + np.sqrt(moments.variance),
        color="red",
    )

    plt.plot(wavelength, moments.mean, color="black")

    plt.xlabel("wavelength [micron]")
    plt.ylabel("normalized intensity [erg cm]")
    plt.show()

    plt.stairs(100 * h / sum(h), edges=locations)
    plt.axvline(locations[imax], color="blue", linestyle="dashed")
    plt.axvline(avg.mean, color="red", linestyle="dashed")
    plt.axvline(med, color="green", linestyle="dashed")
    plt.axvline(avg.mean - 0.5 * np.sqrt(avg.variance), color="blue")
    plt.axvline(avg.mean + 0.5 * np.sqrt(avg.variance), color="blue")

    ax = plt.gca()

    plt.text(
        0.2, 0.70, f"peak: {locations[imax]:0.4f}", color="blue", transform=ax.transAxes
    )

    plt.text(0.2, 0.65, f"<r>: {avg.mean:0.4f}", color="red", transform=ax.transAxes)

    plt.text(
        0.2, 0.60, f"median$_{{r}}$: {med:0.4f}", color="green", transform=ax.transAxes
    )

    plt.text(
        0.2,
        0.55,
        f"$\\sigma_{{\\rm r}}$: {np.sqrt(avg.variance):0.4f}",
        color="yellow",
        transform=ax.transAxes,
    )

    plt.text(
        0.2, 0.50, f"binsize: {locations[1]-locations[0]:0.6f}", transform=ax.transAxes
    )

    plt.text(0.2, 0.45, f"n$_{{\\rm samples}}$: {n}", transform=ax.transAxes)

    plt.text(0.2, 0.40, f"fwhm: {fwhm:0.2f} cm$^{{-1}}$", transform=ax.transAxes)

    plt.text(
        0.2,
        0.35,
        f"E$_{{\\rm in}}$: {Ein / 1.6021765e-12:0.2f} eV",
        transform=ax.transAxes,
    )

    plt.text(
        avg.mean,
        2.5,
        f"{avg.mean:0.4f}",
        transform=ax.transAxes,
        horizontalalignment="center",
    )

    plt.text(
        locations[imax],
        3.5,
        f"{locations[imax]:0.4f}",
        transform=ax.transAxes,
        horizontalalignment="center",
    )

    plt.text(
        med, 4.5, f"{med:0.4}", transform=ax.transAxes, horizontalalignment="center"
    )

    plt.text(
        avg.mean,
        5.6,
        f"{np.sqrt(avg.variance):0.4f}",
        transform=ax.transAxes,
        horizontalalignment="center",
    )

    plt.xlabel("r [Pearson correlation coefficient]")
    plt.ylabel("probability [%]")
    plt.show()

    plt.stairs(100 * hC / sum(hC), edges=locationsC)
    plt.stairs(100 * hDbC / sum(hDbC), edges=locationsDbC, linestyle="dashed")

    plt.axvline(locationsC[imaxC], color="blue")

    plt.axvline(avgC.mean, color="red")
    plt.axvline(avgC.mean - np.sqrt(avgC.variance), color="blue", linestyle="dashed")
    plt.axvline(avgC.mean + np.sqrt(avgC.variance), color="blue", linestyle="dashed")
    plt.plot(
        avgC.mean + np.sqrt(avgC.variance) * np.array([-1, 1]),
        [5.5, 5.5],
        linestyle="dashed",
        color="blue",
    )

    plt.xlabel("<n$_{\\rm carbon}$> [#]")
    plt.ylabel("probability [%]")
    plt.show()

    plt.stairs(100 * hH / sum(hH), edges=locationsH)
    plt.stairs(100 * hDbH / sum(hDbC), edges=locationsDbH, linestyle="dashed")

    plt.axvline(locationsH[imaxH], color="blue")

    plt.axvline(avgH.mean, color="red")
    plt.axvline(avgH.mean - np.sqrt(avgH.variance), color="blue", linestyle="dashed")
    plt.axvline(avgH.mean + np.sqrt(avgH.variance), color="blue", linestyle="dashed")
    plt.plot(
        avgH.mean + np.sqrt(avgH.variance) * np.array([-1, 1]),
        [5.5, 5.5],
        linestyle="dashed",
        color="blue",
    )

    plt.xlabel("<n$_{\\rm hydrogen}$> [#]")
    plt.ylabel("probability [%]")
    plt.show()

    plt.stairs(100 * hCharge / sum(hCharge), edges=locationsCharge)

    plt.axvline(locationsCharge[imaxCharge], color="blue")

    plt.axvline(avgCharge.mean, color="red")
    plt.axvline(
        avgCharge.mean - np.sqrt(avgCharge.variance), color="blue", linestyle="dashed"
    )
    plt.axvline(
        avgCharge.mean + np.sqrt(avgCharge.variance), color="blue", linestyle="dashed"
    )
    plt.plot(
        avgCharge.mean + np.sqrt(avgCharge.variance) * np.array([-1, 1]),
        [5.5, 5.5],
        linestyle="dashed",
        color="blue",
    )

    plt.xlabel("<f$_{\\rm ionized}$> [ratio]")
    plt.ylabel("probability [%]")
    plt.show()
