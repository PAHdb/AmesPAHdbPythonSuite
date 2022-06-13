#!/usr/bin/env python3
"""
test_laboratory.py

Test the laboratory.py module.
"""

import pytest
from pkg_resources import resource_filename

import amespahdbpythonsuite

from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite import laboratory


@pytest.fixture(scope="module")
def pahdb_laboratory():
    xml = "resources/pahdb-experimental_cutdown.xml"
    pahdb = AmesPAHdb(
        filename=resource_filename("amespahdbpythonsuite", xml),
        check=False,
        cache=False,
        update=False,
    )
    return pahdb


class TestLaboratory:
    """
    Test Laboratory class.

    """

    def test_instance(self):
        assert isinstance(
            laboratory.Laboratory(), amespahdbpythonsuite.laboratory.Laboratory
        )

    def test_getset(self, pahdb_laboratory):
        db = pahdb_laboratory
        lab1 = db.getlaboratorybyuid([273])
        l1 = lab1.get()
        assert l1["type"] == "Laboratory"
        lab2 = laboratory.Laboratory()
        lab2.set(l1)
        l2 = lab2.get()
        assert l2["type"] == "Laboratory"
