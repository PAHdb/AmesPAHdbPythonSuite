#!/usr/bin/env python3
"""
test_laboratory.py

Test the laboratory.py module.
"""
import os
import importlib_resources
import matplotlib.pyplot as plt
import pytest

from amespahdbpythonsuite import laboratory
from amespahdbpythonsuite.amespahdb import AmesPAHdb


@pytest.fixture(scope="module")
def test_laboratory():
    xml = (
        importlib_resources.files("amespahdbpythonsuite")
        / "resources/pahdb-experimental_cutdown.xml"
    )
    db = AmesPAHdb(
        filename=xml,
        check=False,
        cache=False,
        update=False,
    )
    return db.getlaboratorybyuid([273])


@pytest.fixture(scope="module")
def test_path(tmp_path_factory):
    d = tmp_path_factory.mktemp("test_observation")
    return f"{d}/result"


class TestLaboratory:
    """
    Test Laboratory class.

    """

    def test_instance(self):
        assert isinstance(laboratory.Laboratory(), laboratory.Laboratory)

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

    def test_write_laboratory(self, test_laboratory, test_path):
        test_laboratory.write(f"{test_path}.tbl")
        assert os.path.exists(f"{test_path}.tbl")
