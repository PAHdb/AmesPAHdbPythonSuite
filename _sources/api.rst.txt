.. sectnum::
   :start: 3

===
API
===

The suite currently consists of 13 IDL object classes and two
helper programs. The main interface to the database is provided by
the *AmesPAHdbIDLSuite*-object, which will parse the dowloaded
XML-file utilizing IDL's AmesPAHXMLParser-object. The available
PAH data is divided into four parts, each handled by its own
object class and providing its own functionality.

1. The Species-object, holds the PAH molecular properties,
2. The Transitions-object, holds the fundamental vibrational transitions,
3. The Laboratory_Spectrum-object holds the raw laboratory spectrum, and
4. The Geometry-objects, holds the geometrical data.

The *AmesPAHdbIDLSuite*-object will return one of these object
classes upon request. In addition, the
*AmesPAHdbIDLSuite*\ \_Transitions-object returns object classes
upon convolving, co-adding and fitting spectra; an
*AmesPAHdbIDLSuite*\ \_Spectrum-,
*AmesPAHdbIDLSuite*\ \_Coadded_Spectrum- and
*AmesPAHdbIDLSuite*\ \_Fitted_Spectrum-object, respectively. The
*AmesPAHdbIDLSuite*\ \_Observation-object can handle some commonly
used astronomical data formats. The
*AmesPAHdbIDLSuite*\ \_Browser-object is a graphical user
interface to the database that can be used to browse its contents.
Lastly, the IDL_Object-object is a class used to inherit from in
IDL versions before 8.0, allowing operator overloaded classes to
still compile on these earlier versions.

The application programming interface (API) and its use is
described in turn.

.. toctree::
   :maxdepth: 1

   sec_03a
   sec_03b
   sec_03c
   sec_03d
   sec_03e
   sec_03f
   sec_03g
   sec_03h
   sec_03i
