"""
example.py


Test the new parser.
"""

from pkg_resources import resource_filename
from ipdb import set_trace as st
from amespahdbpythonsuite import new_xml_parser

from lxml import etree


sample_xml_file = resource_filename(
    'amespahdbpythonsuite', 'resources/pahdb-theoretical_cutdown.xml'
)


test = new_xml_parser.parsey(sample_xml_file)
