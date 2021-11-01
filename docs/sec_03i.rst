.. sectnum::
   :start: 9
   :prefix: 3.

==============================
Spectroscopic Database Fitting
==============================

The first step in fitting astronomical observations is loading
      astronomical observations.

--------------------------------------
Dealing with astronomical observations
--------------------------------------

Astronomical observations can be handled by the
'AmesPAHdbIDLSuite_Observation'-object, which is able to read
text, *ISO*-SWS, and *Spitzer*-IRS-files. A convenience routine is
available to manage units;
AmesPAHdbIDLSuite_CREATE_OBSERVATION_UNITS_S.

.. code:: idl

   observation = OBJ_NEW('AmesPAHdbIDLSuite_Observation', $
                         'myObservationFile', $
                 Units=AmesPAHdbIDLSuite_CREATE_OBSERVATION_UNITS_S())

In addition, an 'AmesPAHdbIDLSuite_Observation'-object can be
created using keyword-initializers.

.. code:: idl

   observation = OBJ_NEW('AmesPAHdbIDLSuite_Observation',
                          X=frequency, $
                          Y=intensity, $
                          ErrY=ystdev, $
                          Units=AmesPAHdbIDLSuite_CREATE_OBSERVATION_UNITS_S())

The 'AmesPAHdbIDLSuite_Observation'-object exposes the observation
and provides the 'Plot', and 'Write'-methods for output. The
'Plot'-method will display the observation and accepts the
'Oplot', and 'Color'-keywords to control overplotting and color,
respectively. Through IDL's keyword inheritance mechanism
additional keywords accepted by IDL's 'PLOT'-procedure can be
passed.

.. code:: idl

   observation->Plot,XRANGE=[2.5,15],/XSTYLE

The 'Write'-method will write the observation to a single text
(.txt) file. Optionally, a prefix can be given that will be
prepended to the filename.

.. code:: idl

   observation->Write,'myPrefix'

**NB** Text-files can have up to five columns, organized as
follows. Column 1: abscissa, Column 2: ordinate, Column 3:
continuum, Column 4: uncertainty in ordinate, and Column 5:
uncertainty in abscissa.

The 'AmesPAHdbIDLSuite_Observation'-object's 'Rebin',
'AbscissaUnitsTo', and 'SetGridRange'-methods can rebin the
observation onto a specified grid or, with the 'Uniform'-Keyword
set, onto a uniform created grid with specified sampling; convert
the units associated with the abscissa; and change the grid range,
respectively.

.. code:: idl

   observation->Rebin,myGrid

   observation->Rebin,5D,/Uniform

   observation->AbsciccaUnitsTo

   observation->SetGridRange,500,2000

----------------
Database fitting
----------------

Spectroscopic database fitting is handled by the
'AmesPAHdbIDLSuite_Spectrum'-object and it either accepts an
'AmesPAHdbIDLSuite_Observation'-object, or a simple array of
ordinates with an optional array of ordinate uncertainties.
Whether ordinate uncertainties are provided or not, the
'AmesPAHdbIDLSuite_Spectrum'-object's 'Fit'-method will perform a
non-negative least-chi-square or non-negative least-square fit and
return an 'AmesPAHdbIDLSuite_Fitted_Spectrum'-object.

.. code:: idl

   fit = spectrum->Fit(intensity, uncertainty)

The 'AmesPAHdbIDLSuite_Fitted_Spectrum'-object exposes the fit and
provides the 'Plot', and 'Write'-methods for output. The
'Plot'-method accepts the 'Residual', 'Size', 'Charge', and
'Composition'-keywords, which selectively display the residual of
the fit, or either the size, charge and compositional breakdown.
Without these keywords the fit itself is displayed.

.. code:: idl

   fit->Plot,/Charge

Optionally, the 'Wavelength', 'Stick', 'Oplot', 'Legend', and
'Color'-keywords can be given to the 'Plot'-method to control the
abscissa, stick representation, overplotting, legend and color,
respectively. Through IDL's keyword inheritance mechanism
additional keywords accepted by IDL's 'PLOT'-procedure can be
passed.

.. code:: idl

   fit->Plot,/Size,/Wavelength,XTITLE=[2.5,15],/XSTYLE

The 'AmesPAHdbIDLSuite_Fitted_Spectrum'-object's 'Write'-method
will write the fit to a single text (.txt) file. Optionally, a
prefix can be given that will be prepended to the filename.

.. code:: idl

   fit->Write,'myPrefix'

The 'AmesPAHdbIDLSuite_Fitted_Spectrum'-object's 'GetClasses', and
'GetBreakdown'-methods return the fit broken down by charge, size,
and composition, where the first provides the spectrum for each
component and the latter its relative contribution.

.. code:: idl

   classes = fit->getClasses()

   breakdown = fit->getBreakdown()

Optionally the 'Small' keyword can be set, which controls the
small cutoff size in number of carbon atoms.

.. code:: idl

   classes = fit->getClasses(Small=20L)

The 'GetBreakdown'-method also accepts the 'Flux'-keyword, which
controls whether the relative breakdown should be reported based
on fitted weight or integrated flux.

.. code:: idl

   breakdown = fit->getBreakdown(/Flux)
