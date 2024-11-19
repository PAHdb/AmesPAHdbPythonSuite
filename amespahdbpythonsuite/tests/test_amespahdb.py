#!/usr/bin/env python3
"""
test_amespahdb.py

Test the amespahdb.py module.
"""
import pytest
from pkg_resources import resource_filename
import amespahdbpythonsuite

from amespahdbpythonsuite.amespahdb import AmesPAHdb


@pytest.fixture(scope="module")
def pahdb_theoretical():
    xml = "resources/pahdb-theoretical_cutdown.xml"
    db = AmesPAHdb(
        filename=resource_filename("amespahdbpythonsuite", xml),
        check=False,
        cache=False,
        update=False,
    )
    return db


@pytest.fixture(scope="module")
def pahdb_laboratory():
    xml = "resources/pahdb-experimental_cutdown.xml"
    db = AmesPAHdb(
        filename=resource_filename("amespahdbpythonsuite", xml),
        check=False,
        cache=False,
        update=False,
    )
    return db


@pytest.fixture(scope="module")
def pahdb_clusters_theoretical():
    xml = "resources/pahdb-clusters-theoretical_cutdown.xml"
    db = AmesPAHdb(
        filename=resource_filename("amespahdbpythonsuite", xml),
        check=False,
        cache=False,
        update=False,
    )
    return db


class TestAmesPAHdb:
    """
    Test AmesPAHdb class.

    """

    def test_instance(self, pahdb_theoretical):
        assert isinstance(pahdb_theoretical, amespahdbpythonsuite.amespahdb.AmesPAHdb)

    def test_update(self):
        xml = "resources/pahdb-theoretical_cutdown.xml"
        path = resource_filename("amespahdbpythonsuite", xml)
        AmesPAHdb(filename=path, check=False, cache=False, update=True)

    def test_type_theoretical(self, pahdb_theoretical):
        assert pahdb_theoretical.gettype() == "theoretical"

    def test_type_laboratory(self, pahdb_laboratory):
        assert pahdb_laboratory.gettype() == "experimental"

    def test_type_clusters_theoretical(self, pahdb_clusters_theoretical):
        assert pahdb_clusters_theoretical.gettype() == "clusters/theoretical"

    def test_version_theoretical(self, pahdb_theoretical):
        assert pahdb_theoretical.getversion() == "3.10"

    def test_version_laboratory(self, pahdb_laboratory):
        assert pahdb_laboratory.getversion() == "2.00"

    def test_version_clusters_theoretical(self, pahdb_clusters_theoretical):
        assert pahdb_clusters_theoretical.getversion() == "1.00"

    def test_env(self):
        # TODO: Turn the sys.exit into exceptions.
        import os

        if "AMESPAHDEFAULTDB" in os.environ:
            del os.environ["AMESPAHDEFAULTDB"]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            AmesPAHdb(check=False, cache=False, update=False)
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 1

    def test_file_not_exist(self):
        # TODO: Turn the sys.exit into exceptions.
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            AmesPAHdb(filename="file_does_not_exist.xml", update=False)
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 2

    def test_cache(self, capsys):
        xml = "resources/pahdb-experimental_cutdown.xml"
        with capsys.disabled():
            AmesPAHdb(
                filename=resource_filename("amespahdbpythonsuite", xml),
                check=False,
                cache=True,
                update=False,
            )
        AmesPAHdb(
            filename=resource_filename("amespahdbpythonsuite", xml),
            check=False,
            cache=True,
            update=False,
        )
        capture = capsys.readouterr()
        assert capture.out.find("RESTORING DATABASE FROM CACHE") >= 0

    def test_checkversion(self, pahdb_theoretical):
        assert not pahdb_theoretical.checkversion("wrong")

    def test_keybyuids(self, pahdb_theoretical):
        uids = [18, 73, 726, 2054, 223]
        res = pahdb_theoretical._AmesPAHdb__getkeybyuids("formula", uids)
        assert sorted(uids) == sorted(list(res.keys()))

    def test_gettransitionsbyuid(self, pahdb_theoretical):
        assert isinstance(
            pahdb_theoretical.gettransitionsbyuid(18),
            amespahdbpythonsuite.transitions.Transitions,
        )

    def test_getlaboratorybyuid(self, pahdb_laboratory):
        assert isinstance(
            pahdb_laboratory.getlaboratorybyuid(273),
            amespahdbpythonsuite.laboratory.Laboratory,
        )

    def test_getlaboratorybyuid_wrong_database(self, pahdb_theoretical):
        assert not pahdb_theoretical.getlaboratorybyuid([273])

    def test_getspeciesbyuid(self, pahdb_theoretical):
        assert isinstance(
            pahdb_theoretical.getspeciesbyuid(18),
            amespahdbpythonsuite.species.Species,
        )

    def test_getgeometrybyuid(self, pahdb_theoretical):
        assert isinstance(
            pahdb_theoretical.getgeometrybyuid(18),
            amespahdbpythonsuite.geometry.Geometry,
        )

    def test_getdatabaseref(self, pahdb_theoretical):
        assert isinstance(pahdb_theoretical.getdatabaseref(), dict)

    def test_search_formula(self, pahdb_theoretical):
        assert pahdb_theoretical.search("c24h12") == [18]

    def test_search_charge(self, pahdb_theoretical):
        assert pahdb_theoretical.search("charge=0") == [18, 73, 726, 2054, 223]

    def test_search_size_and_charge(self, pahdb_theoretical):
        assert pahdb_theoretical.search("c=24 and h=12 and neutral") == [18]

    def test_search_size_ampersand_charge(self, pahdb_theoretical):
        assert pahdb_theoretical.search("c=24 & h=12 & neutral") == [18]

    def test_search_uid_pipe(self, pahdb_theoretical):
        assert pahdb_theoretical.search("uid=18 | uid=73") == [18, 73]

    def test_search_transitions(self, pahdb_theoretical):
        assert pahdb_theoretical.search(
            "c>20 frequency > 3068 and frequency < 3070 intensity > 140 intensity < 141"
        ) == [18]

    def test_search_clusters(self, pahdb_clusters_theoretical):
        assert pahdb_clusters_theoretical.search(
            "monomers=18 type=dimer conformation=step"
        ) == [761058176]
