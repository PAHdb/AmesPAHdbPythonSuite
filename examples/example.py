#!/usr/bin/env python3
"""example.py

Example of using the AmesPAHdbPythonSuite to display the ('stick')
absorption spectrum of coronene (UID=18).

"""

import pkg_resources

from amespahdbpythonsuite.xmlparser import XMLparser

import matplotlib.pyplot as plt

if __name__ == '__main__':

    path = 'resources/pahdb-theoretical_cutdown.xml'

    xml = pkg_resources.resource_filename('amespahdbpythonsuite', path)

    parser = XMLparser(xml)

    parser.verify_schema()

    library = parser.to_pahdb_dict()

    plt.bar([d['frequency'] for d in library['species'][18]['transitions']],
            [d['intensity'] for d in library['species'][18]['transitions']],
            20, color='red', edgecolor="none")

    plt.xlabel('frequency [cm$^{-1}$]')

    plt.ylabel('integrated cross-section [km mol$^{-1}$]')

    plt.show()
