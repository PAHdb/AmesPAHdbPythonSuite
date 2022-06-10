#!/usr/bin/env python3
"""
test_fitted.py

Test the fitted.py module.
"""

import pytest
from os.path import exists
from pkg_resources import resource_filename

from astropy.io import ascii

from amespahdbpythonsuite.amespahdb import AmesPAHdb


@pytest.fixture(scope="module")
def pahdb_theoretical():
    xml = 'resources/pahdb-theoretical_cutdown.xml'
    pahdb = AmesPAHdb(filename=resource_filename('amespahdbpythonsuite', xml),
                      check=False, cache=False, update=False)
    return pahdb


class TestFitted():
    """
    Test Fitted class.

    """

    def test_fit(self, pahdb_theoretical, tmp_path):
        # Read input spectrum.
        spec = resource_filename('amespahdbpythonsuite', 'resources/galaxy_spec.ipac')
        f = ascii.read(spec)
        wave = f['wavelength']
        flux = f['flux']
        sigma = f['sigma']
        waven = [1e4 / x for x in wave]
        # Define output name.
        outputname = spec.split('/')[-1].split('.')[0]
        # Obtain units.
        if f[f.colnames[0]].unit == 'micron':
            xunit = '$\\mu$m'
        else:
            xunit = f[f.colnames[0]].unit
        yunit = f[f.colnames[1]].unit
        units = [xunit, yunit]
        # Read the database.
        pahdb = pahdb_theoretical
        # UIDs test list.
        uids = [18, 73, 726, 2054, 223]
        # Retrieve the transitions from the database for the subset of PAHs.
        transitions = pahdb.gettransitionsbyuid(uids)
        # Set emission model.
        transitions.cascade(6 * 1.603e-12, multiprocessing=False)
        # Shift data 15 wavenumber to the red.
        transitions.shift(-15.0)
        # convolve the transitions into a spectrum.
        spectrum = transitions.convolve(grid=waven, fwhm=15.0, gaussian=True, multiprocessing=False)
        # fit the spectrum.
        fit = spectrum.fit(flux, sigma)
        # Create temporary pytest directory.
        d = tmp_path / "sub"
        d.mkdir()
        out = f'{d}/{outputname}'
        # Create plots.
        fit.plot(wavelength=True, sigma=sigma, outputname=out,
                 ptype='UIDs', ftype='pdf', units=units)
        assert exists(f'{out}_UIDs.pdf')
        fit.plot(wavelength=True, residual=True, sigma=sigma, outputname=out,
                 ptype='residual', ftype='pdf', units=units)
        assert exists(f'{out}_residual.pdf')
        fit.plot(wavelength=True, size=True, sigma=sigma, outputname=out,
                 ptype='size', ftype='pdf', units=units)
        assert exists(f'{out}_size.pdf')
        fit.plot(wavelength=True, charge=True, sigma=sigma, outputname=out,
                 ptype='charge', ftype='pdf', units=units)
        assert exists(f'{out}_charge.pdf')
        fit.plot(wavelength=True, composition=True, sigma=sigma, outputname=out,
                 ptype='composition', ftype='pdf', units=units)
        assert exists(f'{out}_composition.pdf')
        # Save to a temporary file.
        fit.write(out)
        colnames = ['UID', 'formula', 'Nc', 'charge', 'mweight',
                    'n_solo', 'n_duo', 'n_trio', 'n_quartet', 'n_quintet', 'fweight']
        f = ascii.read(f'{out}_results.txt')
        assert f.colnames == colnames
        # Obtain fit breakdown.
        bd = fit.getbreakdown()
        lkeys = ['solo', 'duo', 'trio', 'quartet', 'quintet',
                 'anion', 'neutral', 'cation', 'small', 'large',
                 'nitrogen', 'pure', 'nc', 'err',
                 'e127', 'e112', 'e77', 'e62', 'e33']
        assert list(bd.keys()) == lkeys
