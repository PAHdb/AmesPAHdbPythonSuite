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

The software tools offer three PAH emission models. With increasing complexity they are the 'fixedtemperature', 'calculatedtemperature', and 'cascade' model. The first simply multiplies a blackbody at fixed given temperature with the integrated cross-section of each vibrational transition. The second first calculates the maximum attained temperature from the provided input and subsequently multiplies a blackbody at that fixed temperature with the integrated cross-section of each vibrational transition. The third averages the total emission over the entire cooling cascade (time).

Emission models are handled by the 'transitions'-instance. The 'FixedTemperature'-model simply takes a temperature, in Kelvin, and, in their simplest form, both the 'calculatedtemperature' and 'cascade' models take an energy, in erg.

.. code:: python

   transitions.fixedtemperature(600)

   transitions.calculatedtemperature(4.0 * 1.603e-12)

   transitions.cascade(6 * 1.603e-12)

Both the 'calculatedtemperature' and 'cascade'-methods accept the
'approximate', 'star', 'stellarmodel', and 'isrf'-keywords. With
the 'approximate'-keyword specified, calculations are performed
using the PAH emission model from Bakes et al.
(2001\ `a <http://adsabs.harvard.edu/abs/2001ApJ...556..501B>`__,
`b <http://adsabs.harvard.edu/abs/2001ApJ...560..261B>`__). When
the 'star'-keyword is set, a stellar blackbody at the provided
temperature is used to calculate the average energy absorbed by
each PAH utilizing the PAH absorption cross-sections from Draine &
Li
(`2007 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2007ApJ...657..810D&db_key=AST&link_type=ABSTRACT&high=54888b502c27613>`__).
In case the 'stellarmodel'-keyword is provided as well, the input
is considered to be a full-blown, for example, kurucz stellar
atmosphere model. In the IDL case, the
'AmesPAHdbIDLSuite_CREATE_KURUCZ_STELLARMODEL_S' helper routine is
provided to assist with molding the model data into the proper
input format. Lastly, with the 'isrf'-keyword set, the
interstellar radiation field from Mathis et al.
(`1983 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=1983A%26A...128..212M&db_key=AST&link_type=ABSTRACT&high=54888b502c28367>`__)
is used to calculate the average energy absorbed by each PAH.

.. code:: python

   transitions.calculatedtemperature(17000, star=True)

   transitions.cascade(precision='approximate', field='isrf')

The 'cascade'-method also accepts the 'convolve'-keyword. When set
and combined with either the 'star', optionally with the
'stellarModel'-keyword, or 'isrf'-keyword, will instead of
calculating the average absorbed photon energy for each PAH,
convolve the PAH emission with the entire radiation field.

.. code:: python

   transitions.cascade(17000, star=True, convolve=True)

**NB** This is computationally expensive.

The 'transitions'-instance's 'shift'-method can be
used to redshift the fundamental transitions to simulate
anharmonic effects.

.. code:: python

   transitions.shift(-15)

**NB** Red-shifting the fundamental vibrational transitions should
be done *after* applying one of the three emission models
described above.

-------------
Line profiles
-------------

Line profiles are also handled by the 'transitions'-intance and it provides three profiles Lorentzian, Gaussian and Drude. Convolution with the keyword chosen line profile is achieved through the 'transitions'-instance's 'convolve'-method, which will return the convolved spectrum in the form of a 'spectrum'-instance.

.. code:: python

   spectrum = transitions.convolve(drude=True)

Optionally, the 'convolve'-method accepts the 'fwhm', 'grid',
'npoints', and 'xrange'-keywords, which control the
full-width-at-half-maximum of the selected line profile (in
cm\ :sup:`-1`), convolution onto a specified grid, the number of
resolution elements in the generated spectrum, and the frequency
range (in cm\ :sup:`-1`) of the spectrum.

.. code:: python

   spectrum = transitions.convolve(drude=True, fwhm=20.0, grid=myGrid)

The 'spectrum'-instance exposes convolved spectra and provides the 'plot', and 'write'-methods. The 'plot'-method will display the spectrum of each PAH species in a different color. The 'write'-method will write all spectra to an IPAC table (.tbl). Optionally, a filename can be provided.

.. code:: python

   spectrum.plot()

   spectrum.write('myFile')

Optionally, the 'wavelength', 'stick', 'oplot', 'legend', and
'color'-keywords can be given to the 'plot'-method to control
abscissa, stick representation, overplotting, legend and color,
respectively. 

.. code:: python

   spectrum.plot(wavelength=True)
