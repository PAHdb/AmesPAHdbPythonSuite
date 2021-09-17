.. sectnum::
   :start: 1
   :prefix: 3.

==============
Main Interface
==============

Interaction with the NASA Ames PAH IR Spectroscopic Database is
organized around the 'AmesPAHdbIDLSuite'-object, which is created
as shown below.

.. code:: idl

   pahdb = OBJ_NEW('AmesPAHdbIDLSuite')

Alternatively, the 'AmesPAHdbIDLSuite'-object can be created
through the implicit object creation method as shown below.

.. code:: idl

   pahdb = AmesPAHdbIDLSuite()

In case the system/IDL variable for the default database is not
set, or one wishes to load another version or type of database,
the 'Filename'-keyword can be specified.

.. code:: idl

   pahdb = OBJ_NEW('AmesPAHdbIDLSuite', Filename='/path/to/xml-file')

By default the *AmesPAHdbIDLSuite* will cache the parsed database
XML-file for faster subsequent access. However, this behavior can
be disabled by setting the 'Cache'-keyword to '0'.

.. code:: idl

   pahdb = OBJ_NEW('AmesPAHdbIDLSuite', Cache=0)

When parsing the database XML-file, the *AmesPAHdbIDLSuite* will
try and validate its content against a URL-linked Schema. However,
validation can be disabled by setting the 'Check'-keyword to '0'.
This can be useful when not having an active internet connection.

.. code:: idl

   pahdb = OBJ_NEW('AmesPAHdbIDLSuite', Check=0)

It is possible to combine the different keywords.

.. code:: idl

   pahdb = OBJ_NEW('AmesPAHdbIDLSuite', Filename='path/to/xml-file', $
                                        Cache=0, Check=0)

When finished with the *AmesPAHdbIDLSuite* the object should be
destroyed.

.. code:: idl

   OBJ_DESTROY,pahdb
