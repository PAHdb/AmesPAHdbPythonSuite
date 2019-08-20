=====
Usage
=====

This is the example contained in example.py located in the
examples-directory.

.. code-block:: python

     import pkg_resources
     from amespahdbpythonsuite.xmlparser import XMLparser

     path = 'resources/pahdb-theoretical_cutdown.xml'
     xml = pkg_resources.resource_filename('amespahdbpythonsuite', path)
     parser = XMLparser(xml)
     parser.verify_schema()
     library = parser.to_pahdb_dict()
