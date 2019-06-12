#!/usr/bin/env python3
"""
test_parsers.py

Test whether the XMLparser class behaves as expected.
"""

import pytest
from lxml.etree import XMLSyntaxError
from pkg_resources import resource_filename

import amespahdbpythonsuite
from amespahdbpythonsuite.xmlparser import XMLparser


@pytest.fixture(scope="module")
def real_xml_file():
    xml = 'resources/pahdb-theoretical_cutdown.xml'
    return resource_filename('amespahdbpythonsuite', xml)


@pytest.fixture(scope="module")
def real_xml_file_experimental():
    xml = 'resources/pahdb-experimental_cutdown.xml'
    return resource_filename('amespahdbpythonsuite', xml)


@pytest.fixture(scope="module")
def illformed_xml_file():
    xml = 'resources/pahdb-theoretical_cutdown_illformed.xml'
    return resource_filename('amespahdbpythonsuite', xml)


@pytest.fixture(scope="module")
def not_an_xml_file(tmpdir_factory):
    tmp_file = tmpdir_factory.mktemp("data").join("not_xml.xml")
    with open(tmp_file, 'w') as file:
        file.write('Not an xml file.')
    return tmp_file


class TestXMLparser:
    """Test the XMLparser class."""
    def test_real_xml_file(self, real_xml_file, validate=True):
        parser = XMLparser(real_xml_file, validate=validate)
        parser.verify_schema()
        database, info = parser.to_pahdb_dict()
        assert isinstance(parser, amespahdbpythonsuite.xmlparser.XMLparser)
        assert parser.schema_is_valid is True

        dict_keys = ['filename', 'type', 'date', 'full', 'version',
                     'comment', 'nspecies']
        assert list(info.keys()) == dict_keys
        assert list(database.keys()) == [18]
        xml_str = 'An XML parser from the AmesPAHdbPythonSuite'
        assert str(parser).split('.')[0] == xml_str
        print(parser)

    def test_real_xml_file2(self, real_xml_file, validate=False):
        parser = XMLparser(real_xml_file, validate=validate)
        database, info = parser.to_pahdb_dict()
        assert isinstance(parser, amespahdbpythonsuite.xmlparser.XMLparser)
        assert parser.schema_is_valid is None

        dict_keys = ['filename', 'type', 'date', 'full', 'version',
                     'comment', 'nspecies']
        assert list(info.keys()) == dict_keys
        assert list(database.keys()) == [18]
        print(parser)

    def test_real_xml_file_experimental(self, real_xml_file_experimental):
        parser = XMLparser(real_xml_file_experimental)
        parser.verify_schema()
        database, info = parser.to_pahdb_dict()
        assert isinstance(parser, amespahdbpythonsuite.xmlparser.XMLparser)

        dict_keys = ['filename', 'type', 'date', 'full', 'version',
                     'comment', 'nspecies']
        assert list(info.keys()) == dict_keys
        assert list(database.keys()) == [273]
        assert ('laboratory' in database[273].keys())

    def test_illformed_xml_file(self, illformed_xml_file, validate=True):
        with pytest.raises(XMLSyntaxError):
            parser = XMLparser(illformed_xml_file, validate=validate)
            parser.verify_schema()

    def test_illformed_xml_file2(self, illformed_xml_file, validate=False):
        with pytest.raises(XMLSyntaxError):
            parser = XMLparser(illformed_xml_file, validate=validate)
            parser.verify_schema()
            database, info = parser.to_pahdb_dict()

    def test_not_an_xml_file(self, tmpdir_factory, not_an_xml_file):
        with pytest.raises(TypeError):
            parser = XMLparser(not_an_xml_file)
            parser.verify_schema()
            database, info = parser.to_pahdb_dict()

    with pytest.raises(OSError):
        parser = XMLparser(filename='fake.txt')
        parser.verify_schema()
