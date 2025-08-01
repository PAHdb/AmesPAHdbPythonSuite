#!/usr/bin/env python3

import copy
import io
import re
import sys
from typing import Optional

import astropy.units as u  # type: ignore

from amespahdbpythonsuite import geometry, laboratory, transitions
from amespahdbpythonsuite.amespahdb import AmesPAHdb

message = AmesPAHdb.message


class Species:
    """
    AmesPAHdbPythonSuite species class.
    Contains methods to work with PAH species.

    """

    pahdb = None

    def __init__(self, d: Optional[dict] = None, **keywords) -> None:
        self.set(d, **keywords)

    def set(self, d: Optional[dict] = None, **keywords) -> None:
        """
        Populate properties.

        """
        self.database = keywords.get("database", "")
        self.version = keywords.get("version", "")
        self.data = keywords.get("data", dict())
        self.uids = keywords.get("uids", list())

        if isinstance(d, dict):
            if d.get("type", "") == self.__class__.__name__:
                if "database" not in keywords:
                    self.database = d["database"]
                if "version" not in keywords:
                    self.version = d["version"]
                if "data" not in keywords:
                    self.data = d["data"]
                if "uids" not in keywords:
                    self.uids = d["uids"]

        if "pahdb" in keywords:
            self.pahdb = keywords.get("pahdb")

        if self.pahdb:
            if self.pahdb["database"] != self.database:
                message(
                    f"DATABASE MISMATCH: {self.pahdb['database']} != {self.database}"
                )
                return

            if self.pahdb["version"] != self.version:
                message(f"VERSION MISMATCH: {self.pahdb['version']} != {self.version}")
                return

    def get(self) -> dict:
        """
        Return data dictionary with expected keywords.

        """
        return {
            "type": self.__class__.__name__,
            "database": self.database,
            "version": self.version,
            "data": self.data,
            "uids": self.uids,
        }

    def __repr__(self) -> str:
        """
        Class representation.

        """
        return (
            f"{self.__class__.__name__}({self.uids=},{self.database=},{self.version=})"
        )

    def __str__(self) -> str:
        """
        A description of the instance.
        """

        return f"AmesPAHdbPythonSuite Species instance.\n{self.uids=}"

    def print(self, uid=None, str=False) -> Optional[str]:
        """
        Print species data.

        """
        if uid and uid not in self.data:
            message(f"UID {uid} NOT FOUND")
            return None

        if str:
            sys.stdout = out = io.StringIO()

        if uid:
            print(f"UID       : {uid}")
            for k, v in self.data[uid].items():
                if isinstance(v, list) or isinstance(v, dict):
                    continue
                print(f"{k.upper():<10}: {v}")
            for comment in self.data[uid]["comments"]:
                print(f"COMMENT   : {comment}")
            for reference in self.data[uid]["references"]:
                print(f"REFERENCE : {reference}")
        else:
            for uid in self.uids:
                print("=" * 55)
                print(f"UID       : {uid}")
                for k, v in self.data[uid].items():
                    if isinstance(v, list) or isinstance(v, dict):
                        continue
                    print(f"{k.upper():<10}: {v}")
                for comment in self.data[uid]["comments"]:
                    print(f"COMMENT   : {comment}")
                for reference in self.data[uid]["references"]:
                    print(f"REFERENCE : {reference}")
                print("=" * 55)
        if str:
            sys.stdout = sys.__stdout__
            return out.getvalue()

        return None

    def getuids(self) -> list[int]:
        """
        Return uid list.

        """
        return self.uids

    def intersect(self, uids: list[int]) -> None:
        """
        Updates data to the intersection with provided UIDs.

        Parameters:
            uids : list of integers

        """
        keep = set(self.uids) & set(uids)

        count = len(keep)

        if count == 0:
            message("NO INTERSECTION FOUND")

            return

        message(f"INTERSECTION FOUND: {count}")

        self.uids = list(keep)

        self.data = {key: self.data[key] for key in self.uids}

    def difference(self, uids: list[int]) -> None:
        """
        Updates data to the difference with provided UIDs.

        Parameters:
            uids : list of integers
                List of UIDs.

        """
        keep = set(self.uids) - set(uids)

        count = len(keep)

        if count == 0:
            message("NO DIFFERENCE FOUND")

            return

        message(f"DIFFERENCE FOUND: {keep}")

        self.uids = list(keep)

        self.data = {key: self.data[key] for key in self.uids}

    def transitions(self) -> transitions.Transitions:
        """
        Return transitions instance.
        Calls the :class:`amespahdbpythonsuite.transitions.Transitions` class.

        Returns:
            transitions instance

        """

        return transitions.Transitions(
            database=self.database,
            version=self.version,
            data=self.__getkey("transitions"),
            pahdb=self.pahdb,
            uids=self.uids,
            model={"type": "zerokelvin_m", "temperature": 0.0, "description": ""},
            units={
                "abscissa": {
                    "unit": u.cm**-1,
                    "label": "frequency",
                },
                "ordinate": {"unit": u.km / u.mol, "label": "integrated cross-section"},
            },
        )

    def geometry(self) -> geometry.Geometry:
        """
        Return geometry instance.
        Calls the :class:`amespahdbpythonsuite.geometry.Geometry` class.

        Returns:
            geometry instance

        """

        return geometry.Geometry(
            database=self.database,
            version=self.version,
            data=self.__getkey("geometry"),
            pahdb=self.pahdb,
            uids=self.uids,
        )

    def laboratory(self) -> Optional[laboratory.Laboratory]:
        """
        Return laboratory instance.
        Calls the :class:`amespahdbpythonsuite.laboratory.Laboratory` class.

        Returns:
            laboratory instance

        """

        if self.data["database"] != "experimental":
            message("EXPERIMENTAL DATABASE REQUIRED")
            return None

        return laboratory.Laboratory(
            database=self.database,
            version=self.version,
            data=self.__getkey("laboratory"),
            pahdb=self.pahdb,
            uids=self.uids,
            model={"type": "laboratory_m", "temperature": 0.0, "description": ""},
            units={
                "abscissa": {"unit": 1, "str": "frequency [wavenumber]"},
                "ordinate": {"unit": 2, "str": "absorbance" + "[-log(I/I$_{0})$"},
            },
        )

    def references(self) -> dict:
        """
        Return dict with references.

        Returns:
            dict

        """

        return self.__getkey("references")

    def comments(self) -> dict:
        """
        Return dict with comments.

        Returns:
            dict

        """

        return self.__getkey("comments")

    def __getkey(self, key) -> dict:
        """
        Get a dictionary of species properties
        retrieved by keyword.

        Parameters:
            key : str
                Database keyword.

        Returns:
            Dictionary of retrieved properties with UIDs as keys.

        """

        return copy.deepcopy(
            dict(
                (uid, self.data[uid][key])
                for uid in self.uids
                if uid in self.data.keys()
            )
        )


def formatformula(formula: str) -> str:
    """
    Make the formulae look pretty by embedding LaTeX formatting commands.
    """

    formatted = re.sub(r"([A-Z][a-z]?)([0-9]+)", r"\1$_{\\mathregular{\2}}", formula)

    return re.sub(
        r"((\+)+|(\+[0-9])|(-)+|(-[0-9]))", r"$^{\\mathregular{\1}}", formatted
    )
