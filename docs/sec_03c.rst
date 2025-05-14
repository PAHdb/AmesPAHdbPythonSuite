.. sectnum::
   :start: 3
   :prefix: 3.

=====================================
Working with Molecular PAH Properties
=====================================

The 'species'-instance exposes the molecular PAH properties.

The below code snippets are examples only of the 'AmesPAHdbPythonSuite Species' module. 

.. code:: python

   pahs = pahdb.getspeciesbyuid(pahdb.search("c<=20 neutral n=2 neutral"))

The 'species'-instance's 'print'-method will print
out the associated molecular properties for each PAH species.

.. code:: python

   pahs.print()

Optionally, the 'str'-keyword can be given to the 'print'-method,
which will return the molecular PAH properties for each species as
an array of strings.

.. code:: python

   pahs.print(str=True)
