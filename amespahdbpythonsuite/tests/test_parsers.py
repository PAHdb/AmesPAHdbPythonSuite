#!/usr/bin/env python3
"""
test_parsers.py

Test whether the XMLparser class behaves as expected.
"""

import pandas as pd
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
def illformed_xml_file():
    xml = 'resources/pahdb-theoretical_cutdown_illformed.xml'
    return resource_filename('amespahdbpythonsuite', xml)


@pytest.fixture(scope="module")
def not_an_xml_file(tmpdir_factory):
    tmp_file = tmpdir_factory.mktemp("data").join("not_xml.xml")
    mock_df = pd.DataFrame({'wave': (0, 1), 'flux': (1, 2)})
    mock_df.to_pickle(tmp_file)
    return tmp_file


class TestXMLparser:
    """Test the XMLparser class with a cut-down, real PAHdb XML file."""

    def test_real_xml_file(self, real_xml_file):
        parsey = XMLparser(real_xml_file)
        assert isinstance(parsey, amespahdbpythonsuite.parsers.XMLparser)

        dict_keys = ['filename', 'type', 'date', 'full', 'version',
                     'comment', 'nspecies']
        assert list(parsey.to_dict()[1].keys()) == dict_keys
        assert list(parsey.to_dict()[0].keys()) == [18]

    def test_illformed_xml_file(self, illformed_xml_file):
        with pytest.raises(XMLSyntaxError):
            XMLparser(illformed_xml_file)

    def test_not_an_xml_file(self, tmpdir_factory, not_an_xml_file):
        with pytest.raises(TypeError):
            XMLparser(not_an_xml_file)

    with pytest.raises(OSError):
        XMLparser('fake.txt')
