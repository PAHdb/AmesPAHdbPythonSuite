#!/usr/bin/env python3

from amespahdbpythonsuite import laboratory
from amespahdbpythonsuite import geometry
from amespahdbpythonsuite import transitions
import copy
import re

from typing import Union

from amespahdbpythonsuite.amespahdb import AmesPAHdb

message = AmesPAHdb.message


class Species:
    """
    AmesPAHdbPythonSuite species class.
    Contains methods to work with PAH species.

    """

    def __init__(self, d=None, **keywords):
        self.type = ""
        self.version = ""
        self.data = dict()
        self.pahdb = None
        self.uids = []
        self.set(d, **keywords)

    def set(self, d, **keywords) -> None:
        """
        Populate properties.

        """
        if d:
            if d.get("type", "") == self.__class__.__name__:
                if not keywords.get("type"):
                    self.type = d["database"]
                if not keywords.get("version"):
                    self.version = d["version"]
                if not keywords.get("data"):
                    self.data = d["data"]
                if not keywords.get("uids"):
                    self.uids = d["uids"]

        if keywords.get("type"):
            self.type = keywords.get("type")
        if keywords.get("version"):
            self.version = keywords.get("version")
        if keywords.get("data"):
            self.data = keywords.get("data")
        if keywords.get("pahdb"):
            self.pahdb = keywords.get("pahdb")
        if keywords.get("uids"):
            self.uids = keywords.get("uids")

        if self.pahdb:
            if self.pahdb["database"] != self.type:

                message(f'DATABASE MISMATCH: {self.pahdb["database"]} != {self.type}')
                return

            if self.pahdb["version"] != self.version:

                message(f'VERSION MISMATCH: {self.pahdb["version"]} != {self.version}')
                return

    def get(self) -> dict:
        """
        Return data dictionary with expected keywords.

        """
        return {
            "type": self.__class__.__name__,
            "database": self.type,
            "version": self.version,
            "data": self.data,
            "uids": self.uids,
        }

    def getuids(self) -> list:
        """
        Return uid list.

        """
        return self.uids

    def intersect(self, uids: list) -> None:
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

        self.uids = keep

        self.data = {key: self.data[key] for key in self.uids}

    def difference(self, uids: list) -> None:
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

        self.uids = keep

        self.data = {key: self.data[key] for key in self.uids}

    def transitions(self) -> transitions.Transitions:
        """
        Return transitions object.
        Calls the :class:`amespahdbpythonsuite.transitions.Transitions` class.

        Returns:
            transitions object

        """

        return transitions.Transitions(
            type=self.type,
            version=self.version,
            data=self.__getkey("transitions"),
            pahdb=self.pahdb,
            uids=self.uids,
            model={"type": "zerokelvin_m", "temperature": 0.0, "description": ""},
            units={
                "abscissa": {"unit": 1, "str": "frequency [wavenumber]"},
                "ordinate": {"unit": 2, "str": "integrated cross-section" + "[km/mol]"},
            },
        )

    def geometry(self) -> geometry.Geometry:
        """
        Return geometry object.
        Calls the :class:`amespahdbpythonsuite.geometry.Geometry` class.

        Returns:
            geometry object

        """

        return geometry.Geometry(
            type=self.type,
            version=self.version,
            data=self.__getkey("geometry"),
            pahdb=self.pahdb,
            uids=self.uids,
        )

    def laboratory(self) -> Union[laboratory.Laboratory, None]:
        """
        Return laboratory object.
        Calls the :class:`amespahdbpythonsuite.laboratory.Laboratory` class.

        Returns:
            laboratory object

        """

        if self.__data["database"] != "experimental":

            message("EXPERIMENTAL DATABASE REQUIRED")

            return None

        return laboratory.Laboratory(
            type=self.type,
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
