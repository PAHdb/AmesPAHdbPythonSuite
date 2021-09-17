.. sectnum::
   :start: 6
   :prefix: 3.

===================================
Working with Raw Laboratory Spectra
===================================

The 'AmesPAHdbIDLSuite_Laboratory'-object exposes available raw
laboratory spectra when an experimental database XML-file is
loaded.

.. code:: idl

   laboratory = pahdb->getLaboratoryByUID( $
                pahdb->Search("c<=20 neutral n=2 neutral"))

The 'AmesPAHdbIDLSuite_Laboratory-object's 'Plot'-method will
display the raw laboratory spectra. The spectrum of each PAH
species will be presented in a different color.

.. code:: idl

   laboratory->Plot

Optionally, the 'Wavelength', 'Oplot', and 'Color'-keywords can be
given to the 'Plot'-method to control the abscissa, overplotting,
and color, respectively. Through IDL's keyword inheritance
mechanism additional keywords accepted by IDL's 'PLOT'-procedure
can be passed.

.. code:: idl

   laboratory->Plot,/Wavelength,Color=5,XRANGE=[2.5,15],/XSTYLE

The 'AmesPAHdbIDLSuite_Laboratory-object's 'Write'-method will
write the raw laboratory spectrum to file. The spectrum of each
PAH species will be written to a separate text (.txt) file, where
each filename will have the PAH species UID embedded. Optionally,
a prefix can be given that will be prepended to the filename.

.. code:: idl

   laboratory->Write,'myPrefix'
