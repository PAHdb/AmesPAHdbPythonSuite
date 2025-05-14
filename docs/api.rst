.. sectnum::
   :start: 3

===
API
===

The main interface to the database is provided by the *AmesPAHdb* instance,
which will parse the dowloaded XML-file utilizing a *XMLParser* instance. The
available PAH data is divided into four parts, each handled by its own class
and providing its own functionality.

1. The *Species* class, which deals with the PAH molecular properties,
2. The *Transitions* class, which deals with the fundamental vibrational transitions,
3. The *Laboratory* class, which deals with the raw laboratory spectrum, and
4. The *Geometry* class, which deals with the geometrical data.

The *AmesPAHdb* instance will return one of these class instances upon request.
In addition, the *Transitions* instance returns class instances upon
convolving, co-adding and fitting spectra; a *Spectrum*, *Coadded* and *Fitted*
or *MCFitted* instance, respectively. The *Observation* class can handle some
commonly used astronomical data formats.

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
