.. sectnum::
   :start: 8
   :prefix: 3.

====================
Creating PAH Spectra
====================

The (theoretical) database XML-files provide fundamental
vibrational transitions at 0-Kelvin. To compare these with
observations, they need to be transformed into a spectral density,
i.e., spectra. When not dealing with absorption at 0-Kelvin, an
emission model is required as well.

---------------
Emission models
---------------

The *AmesPAHdbIDLSuite* offers three PAH emission models. With
increasing complexity they are the 'FixedTemperature',
'CalculatedTemperature', and 'Cascade' model. The first simply
multiplies a blackbody at fixed given temperature with the
integrated cross-section of each vibrational transition. The
second first calculates the maximum attained temperature from the
provided input and subsequently multiplies a blackbody at that
fixed temperature with the integrated cross-section of each
vibrational transition. The third averages the total emission over
the entire cooling cascade (time).

Emission models are handled by the
'AmesPAHdbIDLSuite_Transitions'-object. The
'FixedTemperature'-model simply takes a temperature, in Kelvin,
and, in their simplest form, both the 'CalculatedTemperature' and
'Cascade' models take an energy, in erg.

.. code:: idl

   transitions->FixedTemperature,600D

   transitions->CalculatedTemperature,6D*1.603D-12; 6 eV

   transitions->Cascade,6D*1.603D-12; 6 eV

Both the 'CalculatedTemperature' and 'Cascade'-methods accept the
'Approximate', 'Star', 'StellarModel', and 'ISRF'-keywords. With
the 'Approximate'-keyword specified, calculations are performed
using the PAH emission model from Bakes et al.
(2001\ `a <http://adsabs.harvard.edu/abs/2001ApJ...556..501B>`__,
`b <http://adsabs.harvard.edu/abs/2001ApJ...560..261B>`__). When
the 'Star'-keyword is set, a stellar blackbody at the provided
temperature is used to calculate the average energy absorbed by
each PAH utilizing the PAH absorption cross-sections from Draine &
Li
(`2007 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2007ApJ...657..810D&db_key=AST&link_type=ABSTRACT&high=54888b502c27613>`__).
In case the 'StellarModel'-keyword is provided as well, the input
is considered to be a full-blown, for example, Kurucz stellar
atmosphere model. The
'AmesPAHdbIDLSuite_CREATE_KURUCZ_STELLARMODEL_S' helper routine is
provided to assist with molding the model data into the proper
input format. Lastly, with the 'ISRF'-keyword set, the
interstellar radiation field from Mathis et al.
(`1983 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=1983A%26A...128..212M&db_key=AST&link_type=ABSTRACT&high=54888b502c28367>`__)
is used to calculate the average energy absorbed by each PAH.

.. code:: idl

   transitions->CalculatedTemperature,17D3,/Star

   transitions->Cascade,/Approximate,/ISRF

   FTAB_EXT,'ckp00_17000.fits',[1,10],angstroms,flam,EXT=1

   transitions->Cascade, $
                AmesPAHdbIDLSuite_CREATE_KURUCZ_STELLARMODEL_S(angstroms, $
                                                               flam), $
                /Star, $
                /StellarModel

The 'Cascade'-method also accepts the 'Convolve'-keyword. When set
and combined with either the 'Star', optionally with the
'StellarModel'-keyword, or 'ISRF'-keyword, will instead of
calculating the average absorbed photon energy for each PAH,
convolve the PAH emission with the entire radiation field.

.. code:: idl

   transitions->Cascade,17D3,/Star,/Convolve

**NB** This is computationally expensive.

The 'AmesPAHdbIDLSuite_Transitions'-object's 'Shift'-method can be
used to redshift the fundamental transitions to simulate
anharmonic effects.

.. code:: idl

   transitions->Shift,-15D

**NB** Red-shifting the fundamental vibrational transitions should
be done *after* applying one of the three emission models
described above.

-------------
Line profiles
-------------

Line profiles are handled by the
'AmesPAHdbIDLSuite_Transitions'-object and it provides three
profiles; Lorentzian, Gaussian and Drude. Convolution with the
keyword chosen line profile is achieved through the
'AmesPAHdbIDLSuite_Transitions'-object's 'Convolve'-method, which
will return the convolved spectrum in the form of an
'AmesPAHdbIDLSuite_Spectrum'-object.

.. code:: idl

   spectrum = transitions->Convolve(/Drude)

Optionally, the 'Convolve'-method accepts the 'FWHM', 'Grid',
'NPoints', and 'XRange'-keywords, which control the
full-width-at-half-maximum of the selected line profile (in
cm\ :sup:`-1`), convolution onto a specified grid, the number of
resolution elements in the generated spectrum, and the frequency
range (in cm\ :sup:`-1`) of the spectrum.

.. code:: idl

   spectrum = transitions->Convolve(/Drude, FWHM=20D, Grid=myGrid)

The 'AmesPAHdbIDLSuite_Spectrum'-object exposes the convolved
spectra and provides the 'Plot', and 'Write'-methods. The
'Plot'-method will display the spectrum of each PAH species in a
different color. The 'Write'-method will write all spectra to a
single text (.txt) file. Optionally, a prefix can be given that
will be prepended to the filename.

.. code:: idl

   spectrum->Plot

   spectrum->Write,'myPrefix''

 Optionally, the 'Wavelength', 'Stick', 'Oplot', 'Legend', and
 'Color'-keywords can be given to the 'Plot'-method to control
 abscissa, stick representation, overplotting, legend and color,
 respectively. Through IDL's keyword inheritance mechanism
 additional keywords accepted by IDL's 'PLOT'-procedure can be
 passed.

.. code:: idl

   spectrum->Plot,/Wavelength,XRANGE=[2.5,15],/XSTYLE
