.. sectnum::
   :start: 2
   :prefix: 3.

================================
Accessing the molecular PAH data
================================

The molecular PAH data is accessed through their associated unique
identifier (UID). However, a search-interface is available that
allows retrieval of these UIDs through simplified queries. The
syntax is the same as used on the NASA Ames PAH IR Spectroscopic
Database
`website. <https://www.astrochemistry.org/pahdb/theoretical/3.20/default/view>`__

.. code:: python

   uids = pahdb.search("c<=20 neutral n=2 neutral")

The molecular PAH data contained in the database XML-files
consists of four main parts:

1. Molecular properties, .e.g, molecular weight, zero point
   energy, total energy, etc.
2. Fundamental vibrational transitions
3. Molecular geometric data
4. Raw laboratory spectra for a subset of molecules in the
   experimental database

.. code:: python

   pahdb = AmesPAHdb()
   
   uids = pahdb.search("c<=20 neutral n=2 neutral")

   pahs = pahdb.getspeciesbyuid(uids)

   transitions = pahdb.gettransitionsbyuid(uids)

   geometry = pahdb.getgeometrybyuid(uids)

   laboratory = pahdb.getlaboratorybyuid(uids)

Alternatively, one can access these components through the 'species'-instance.

.. code:: python

   pahdb = AmesPAHdb()

   pahs = pahdb.getspeciesbyuid(pahdb.search("c<=20 neutral n=2 neutral"))

   transitions = pahs.transitions()

   geometry = pahs.geometry()

   laboratory = pahs.laboratory()
