#!/usr/bin/env python3

from typing import Optional

from amespahdbpythonsuite.amespahdb import AmesPAHdb

message = AmesPAHdb.message


class Data(object):
    """
    AmesPAHdbPythonSuite data class

    """

    def __init__(self, d: Optional[dict] = None, **keywords) -> None:
        """
        Initialize Data class.

        """
        self.__set(d, **keywords)

    def set(self, d: Optional[dict] = None, **keywords) -> None:
        """
        Populate data dictionary.

        """
        self.__set(d, **keywords)

    def __set(self, d: Optional[dict] = None, **keywords):
        """
        Populate data dictionary helper.

        """
        if isinstance(d, dict):
            # Check if expected keywords are present in provided dictionary,
            # otherwise assign them to instance variables.
            if d.get("type", "") == self.__class__.__name__:
                if "type" not in keywords:
                    self.type = d["database"]
                if "version" not in keywords:
                    self.version = d["version"]
                if "data" not in keywords:
                    self.data = d["data"]
                if "uids" not in keywords:
                    self.uids = d["uids"]
                if "model" not in keywords:
                    self.model = d["model"]
                if "units" not in keywords:
                    self.units = d["units"]

        self.type = keywords.get("type", "")
        self.version = keywords.get("version", "")
        self.data = keywords.get("data", dict())
        self.pahdb = keywords.get("pahdb", None)
        self.uids = keywords.get("uids", list())
        self.model = keywords.get("model", dict())
        self.units = keywords.get("units", dict())

        # Check for database and versioning mismatch between provided dictionary and parsed database.
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
            "model": self.model,
            "units": self.units,
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

        self.uids = list(keep)

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

        self.uids = list(keep)

        self.data = {key: self.data[key] for key in self.uids}
