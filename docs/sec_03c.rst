.. sectnum::
   :start: 3
   :prefix: 3.

=====================================
Working with Molecular PAH Properties
=====================================

The 'AmesPAHdbIDLSuite_Species'-object exposes the molecular PAH
properties.

.. code:: idl

   pahs = pahdb->getSpeciesByUID( $
          pahdb->Search("c<=20 neutral n=2 neutral"))

The 'AmesPAHdbIDLSuite_Species'-object's 'Print'-method will print
out the associated molecular properties for each PAH species.

.. code:: idl

   pahs->Print

Optionally, the 'Str'-keyword can be given to the 'Print'-method,
which will return the molecular PAH properties for each species as
an array of strings.

.. code:: idl

   pahs->Print,Str=Str

