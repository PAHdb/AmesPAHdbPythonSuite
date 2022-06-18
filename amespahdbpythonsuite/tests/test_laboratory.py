#!/usr/bin/env python3
"""
test_laboratory.py

Test the laboratory.py module.
"""

import pytest
from pkg_resources import resource_filename
import matplotlib.pyplot as plt

from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite import laboratory


@pytest.fixture(scope="module")
def test_laboratory():
    xml = "resources/pahdb-experimental_cutdown.xml"
    db = AmesPAHdb(
        filename=resource_filename("amespahdbpythonsuite", xml),
        check=False,
        cache=False,
        update=False,
    )
    return db.getlaboratorybyuid([273])


class TestLaboratory:
    """
    Test Laboratory class.

    """

    def test_instance(self):
        assert isinstance(
            laboratory.Laboratory(), laboratory.Laboratory
        )

    def test_plot(self, monkeypatch, test_laboratory):
        monkeypatch.setattr(plt, "show", lambda: None)
        test_laboratory.plot(show=True)

    def test_getset(self, test_laboratory):
        l1 = test_laboratory.get()
        assert l1["type"] == "Laboratory"
        lab2 = laboratory.Laboratory()
        lab2.set(l1)
        l2 = lab2.get()
        assert l2["type"] == "Laboratory"
