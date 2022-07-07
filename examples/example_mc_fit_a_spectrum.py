#!/usr/bin/env python3
"""example2.py

Example of using the AmesPAHdbPythonSuite to display the ('stick')
absorption spectrum of coronene (UID=18).

For more examples visit the PAHdb cookbook website:
https://pahdb.github.io/cookbook/

"""

from pkg_resources import resource_filename

from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite import observation

if __name__ == '__main__':

    obs = observation.Observation(
        resource_filename("amespahdbpythonsuite", "resources/galaxy_spec.ipac")
    )
    # Read the database.
    xml = 'resources/pahdb-theoretical_cutdown.xml'
    pahdb = AmesPAHdb(filename=resource_filename('amespahdbpythonsuite', xml),
                      check=False, cache=False)

    uids = pahdb.search("magnesium=0 oxygen=0 iron=0 silicium=0 chx=0 ch2=0 c>20 h>0")

    # Put back fullerenes, which have h=0
    fullerenes = [717, 720, 723, 735, 736, 737]

    uids = uids + fullerenes

    # Retrieve the transitions from the database for the subset of PAHs.
    transitions = pahdb.gettransitionsbyuid([uids])

    # Calculate the emission spectrum at the temperature reached
    # after absorbing a 6 eV (CGS units) photon.
    transitions.cascade(6 * 1.603e-12, multiprocessing=False)

    # Shift data 15 wavenumber to the red
    transitions.shift(-15.0)

    # Convolve the bands with a Gaussian with FWHM of 15 /cm.
    spectrum = transitions.convolve(grid=1e4 / obs.getgrid(), fwhm=15.0, gaussian=True, multiprocessing=False)

    # Fit the spectrum using Monte Carlo approach.
    mcfit = spectrum.mcfit(obs, nsamples=10)

    # Get average/min/max/std statistics for the breakdown components and uncertainty, and save to file.
    mcfit.stats(mcfit.mcresults)

    # Get average and std spectra and dump them into pickle.
    mccomponents = mcfit.averagespec(mcfit.mcfits, mcfit.mcclasses, dump=True)

    # Create plots.
    mcfit.plot(mccomponents, wavelength=True, sigma=obs.spectrum.uncertainty.array,
               ptype='charge', ftype='pdf')

    mcfit.plot(mccomponents, wavelength=True, sigma=obs.spectrum.uncertainty.array,
               ptype='size', ftype='pdf')

    mcfit.plot(mccomponents, wavelength=True, sigma=obs.spectrum.uncertainty.array,
               ptype='composition', ftype='pdf')
