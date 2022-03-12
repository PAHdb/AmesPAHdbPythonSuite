.. sectnum::
   :start: 9
   :prefix: 3.

==============================
Spectroscopic Database Fitting
==============================

--------------------------------------
Dealing with astronomical observations
--------------------------------------

Astronomical observations can be handled by the 'observation'-instance, which is able to read text, ISO-SWS, and Spitzer-IRS-files (*work-in-progress*).

.. code:: python

   observation = observation('myObservationFile')

**NB** Text-files can have up to five columns, organized as follows:

Column 1: abscissa

Column 2: ordinate

Column 3: continuum

Column 4: uncertainty in ordinate

Column 5: uncertainty in abscissa

.. code:: idl

   observation = observation(x=frequency, \
                             y=intensity, \
                             yerr=ystdev)

The 'observation'-instance exposes the observation and provides the 'plot', and 'write'-methods for output. The 'plot'- method will display the observation and accepts the 'oplot', and 'color'-keywords to control overplotting and color, respectively.

.. code:: python

   observation.plot()

The 'Write'-method will write the observation to an IPAC table (.tbl) file. Optionally, a filename can be provided.

.. code:: python

   observation.write('myFile')


The 'observations'-instance's 'rebin',
'abscissaunitsto', and 'setgridrange'-methods can rebin the
observation onto a specified grid or, with the 'uniform'-keyword
set, onto a uniform created grid with specified sampling; convert
the units associated with the abscissa; and change the grid range,
respectively.

.. code:: python

   observation.rebin(myGrid)

   observation.rebin(5.0, uniform=True)

   observation.AbsciccaUnitsTo

   observation.setgridrange(500, 2000)

----------------
Database fitting
----------------

Spectroscopic database fitting is handled by the 'spectrum'-instance and it either accepts an 'observation'-instance, or simply an array of ordinates with an optional array of ordinate uncertainties. Whether ordinate uncertainties are provided or not, the 'spectrum'- instance's 'fit'-method will perform a non-negative least-chi-square or non-negative least-square fit and return an 'fit'-instance.

.. code:: python

   # Using an 'observation'-instance
   fit = spectrum.fit(observation)

   # Using an array of ordinate values
   fit = spectrum.fit(intensity)

   # Using an array of ordinate uncertainty values
   fit = spectrum.fit(intensity, uncertainty)

The 'fit'-instance exposes the fit and provides the 'plot', and 'Write'-methods for output. The 'plot'-method accepts the 'residual', 'size', 'charge', and 'composition'-keywords, which selectively display the residual of the fit, or either the size, charge and compositional breakdown. Without these keywords the fit itself is displayed.

.. code:: python

   fit.plot(charge=True)

Optionally, the 'wavelength', 'stick', 'oplot', 'legend', and
'color'-keywords can be given to the 'plot'-method to control the
abscissa, stick representation, overplotting, legend and color,
respectively.

.. code:: python

   fit.plot(size=True, wavelength=True)

The 'fit’-instance’s ‘Write’-method will write the fit to an IPAC table (.tbl) file. Optionally, a filename can be provided.

.. code:: python

   fit.write('myFile')

The 'fit'-instance's 'getclasses', and 'getbreakdown'-methods return the fit broken down by charge, size, and composition, where the first provides the spectrum for each component and the latter its relative contribution.

.. code:: python

   classes = fit.getclasses()

   breakdown = fit.getbreakdown()

Optionally the 'small' keyword can be set, which controls the
small cutoff size in number of carbon atoms.

.. code:: python

   classes = fit.getclasses(small=20)

The 'getbreakdown'-method also accepts the 'flux'-keyword, which controls whether the relative breakdown should be reported based on fitted weight or integrated flux.

.. code:: python

   breakdown = fit.getbreakdown(flux=True)
