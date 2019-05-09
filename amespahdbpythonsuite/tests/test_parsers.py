#!/usr/bin/env python3
"""
test_parsers.py

Test whether the XMLparser class behaves as expected.
"""

import pytest
from lxml.etree import XMLSyntaxError
from pkg_resources import resource_filename

import amespahdbpythonsuite
from amespahdbpythonsuite.parsers import XMLparser


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
    def test_real_xml_file(self, real_xml_file, check_schema=True):
        parsed = XMLparser(real_xml_file)
        assert isinstance(parsed, amespahdbpythonsuite.parsers.XMLparser)

        dict_keys = ['filename', 'type', 'date', 'full', 'version',
                     'comment', 'nspecies']
        assert list(parsed.to_dicts()[1].keys()) == dict_keys
        assert list(parsed.to_dicts()[0].keys()) == [18]
        xml_str = 'An XML parser from the Ames PAHdb Python Suite'
        assert str(parsed).split('.')[0] == xml_str
        print(parsed)

    def test_real_xml_file2(self, real_xml_file, check_schema=False):
        parsed = XMLparser(real_xml_file)
        assert isinstance(parsed, amespahdbpythonsuite.parsers.XMLparser)

        dict_keys = ['filename', 'type', 'date', 'full', 'version',
                     'comment', 'nspecies']
        assert list(parsed.to_dicts()[1].keys()) == dict_keys
        assert list(parsed.to_dicts()[0].keys()) == [18]
        xml_str = 'An XML parser from the Ames PAHdb Python Suite'
        assert str(parsed).split('.')[0] == xml_str

    def test_real_xml_file_experimental(self, real_xml_file_experimental):
        parsed = XMLparser(real_xml_file_experimental)
        assert isinstance(parsed, amespahdbpythonsuite.parsers.XMLparser)

        dict_keys = ['filename', 'type', 'date', 'full', 'version',
                     'comment', 'nspecies']
        assert list(parsed.to_dicts()[1].keys()) == dict_keys
        assert list(parsed.to_dicts()[0].keys()) == [273]
        assert ('laboratory' in parsed.to_dicts()[0][273].keys())

    def test_illformed_xml_file(self, illformed_xml_file, check_schema=True):
        with pytest.raises(XMLSyntaxError):
            XMLparser(illformed_xml_file)

    def test_illformed_xml_file2(self, illformed_xml_file, check_schema=False):
        with pytest.raises(XMLSyntaxError):
            XMLparser(illformed_xml_file)

    def test_not_an_xml_file(self, tmpdir_factory, not_an_xml_file):
        with pytest.raises(TypeError):
            XMLparser(not_an_xml_file)

    with pytest.raises(OSError):
        XMLparser('fake.txt')
