.. sectnum::
   :start: 1
   :prefix: 3.

==============
Main Interface
==============

Interaction with the NASA Ames PAH IR Spectroscopic Database is
organized around the 'AmesPAHdb'-object, which is created
as shown below.

.. code:: python

   from amespahdbpythonsuite.amespahdb import AmesPAHdb
   pahdb = AmesPAHdb()

In case the system variable for the default database is not
set, or one wishes to load another version or type of database,
the 'filename'-keyword can be specified.

.. code:: python

   pahdb = AmesPAHdb(filename='/path/to/xml-file')

By default the parsed database XML-file will be cached for faster subsequent access. However, this behavior can be disabled by setting the 'cache'-keyword to false.

.. code:: python

   pahdb = AmesPAHdb(cache=False)

When parsing a database XML-file, the software will validate its content against a URL-linked Schema. However, validation can be disabled by setting the 'check'-keyword to false. This can be useful when not having an active internet connection.

.. code:: python

   pahdb = AmesPAHdb(check=False)

It is possible to combine the different keywords.

.. code:: python

   pahdb = AmesPAHdb(filename='/path/to/xml-file', cache=False, check=False)

Lastly, when finished with the ‘pahdb’-instance it should be destroyed when no garbage collection is available.

.. code:: python

   del pahdb
