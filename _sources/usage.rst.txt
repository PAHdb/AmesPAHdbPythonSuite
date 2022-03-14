.. sectnum::
   :start: 2

=====
Usage
=====

This is the example contained in example.py located in the
examples-directory and shows how the AmesPAHdbPythonSuite is used to
display the ('stick') absorption spectrum of coronene (UID=18).

.. code-block:: python

        from pkg_resources import resource_filename
        from amespahdbpythonsuite.amespahdb import AmesPAHdb

        # Read the database.
        xml = 'resources/pahdb-theoretical_cutdown.xml'
        pahdb = AmesPAHdb(filename=resource_filename('amespahdbpythonsuite', xml), check=False, cache=False)

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
