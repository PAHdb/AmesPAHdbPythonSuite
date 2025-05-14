.. sectnum::
   :start: 6
   :prefix: 3.

===================================
Working with Raw Laboratory Spectra
===================================

The 'laboratory'-instance exposes available raw
laboratory spectra when an experimental database XML-file is
loaded.

.. code:: python

   uids = pahdb.search("c<=20 neutral n=2 neutral")

   laboratory = pahdb.getlaboratorybyuid(uids)


The 'laboratory'-instance's 'plot'-method will
display the raw laboratory spectra. The spectrum of each PAH
species will be presented in a different color.

.. code:: python

   laboratory.plot()

Optionally, the 'wavelength', 'oplot', and 'color'-keywords can be
given to the 'plot'-method to control the abscissa, overplotting,
and color, respectively.

.. code:: python

   laboratory.plot(wavelength=True,color='blue')

The 'laboratory'-instance's 'write'-method will
write the raw laboratory spectrum to file. The spectrum of each
PAH species will be written to a separate text (.txt) file, where
each filename will have the PAH species UID embedded. Optionally,
a prefix can be given that will be prepended to the filename.

.. code:: python

   laboratory.write('myFile')
