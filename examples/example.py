#!/usr/bin/env python3
"""example.py

Example of using the AmesPAHdbPythonSuite.

"""

import pkg_resources

from amespahdbpythonsuite.xmlparser import XMLparser

if __name__ == '__main__':

    path = 'resources/pahdb-theoretical_cutdown.xml'

    xml = pkg_resources.resource_filename('amespahdbpythonsuite', path)

    parser = XMLparser(xml)

    parser.verify_schema()

    library = parser.to_pahdb_dict()
