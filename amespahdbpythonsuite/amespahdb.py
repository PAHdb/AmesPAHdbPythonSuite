#!/usr/bin/env python3

from __future__ import annotations

import copy
import hashlib
import json
import os
import pickle
import random
import re
import sys
import tempfile
import time
import urllib.request
from datetime import timedelta
from typing import TYPE_CHECKING, Optional, Union
from urllib.error import HTTPError

import astropy.units as u  # type: ignore
import numpy as np
from packaging.version import Version

import amespahdbpythonsuite as suite
from amespahdbpythonsuite.xmlparser import XMLparser

if TYPE_CHECKING:
    from amespahdbpythonsuite.geometry import Geometry
    from amespahdbpythonsuite.laboratory import Laboratory
    from amespahdbpythonsuite.species import Species
    from amespahdbpythonsuite.transitions import Transitions


class AmesPAHdb:
    """
    AmesPAHdbPythonSuite main class.
    Contains methods to parse the database, perform search based on query,
    and retrieve UIDs.
    Calls classes:
    :class:`amespahdbpythonsuite.transitions.Transitions`,
    :class:`amespahdbpythonsuite.laboratory.Laboratory`,
    :class:`amespahdbpythonsuite.species.Species`,
    :class:`amespahdbpythonsuite.geometry.Geometry`,
    to retrieve the respective instances.

    """

    def __init__(self, **keywords) -> None:
        """
        Initialize amespahdbpythonsuite class. Prints basic PAHdb info
        and calls the :meth:`amespahdbpythonsuite.amespahdb.parsefile`
        method to parse or restore the database.

        """

        intro = [
            "AmesPAHdbPythonSuite\n",
            "by\n",
            "Dr. Christiaan Boersma\n",
            "and\n",
            "Dr. Alexandros Maragkoudakis\n",
            "Dr. Matthew J. Shannanon\n",
            "Dr. Joseph E. Roser\n",
        ]
        self.message(intro)

        self.message(f"SUITE VERSION: {suite.__version__}")

        if keywords.get("update", False) or (
            keywords.get("update", True) and random.randint(0, 4) == 4
        ):
            self.message("CHECKING FOR UPDATE")
            github = "http://api.github.com/repos/pahdb/amespahdbpythonsuite/tags"
            try:
                with urllib.request.urlopen(github) as url:
                    data = json.load(url)
                    versions = [Version(tag["name"]) for tag in data]
                    versions.sort()
                    if Version(suite.__version__) < versions[-1]:
                        update = versions[-1].public
                        self.message(f"V{update} UPDATE AVAILABLE")
                    else:
                        self.message("NO UPDATE AVAILABLE")
            except HTTPError:
                self.message("FAILED TO CHECK FOR UPDATE")

        self.message("WEBSITE: WWW.ASTROCHEM.ORG/PAHDB/")
        self.message("CONTACT: CHRISTIAAN.BOERSMA@NASA.GOV")

        filename = keywords.get("filename")
        if not filename:
            filename = os.environ.get("AMESPAHDEFAULTDB")
            if not filename:
                msg = [
                    "DATABASE NOT FOUND:",
                    "SET SYSTEM AMESPAHDEFAULTDB",
                    "ENVIRONMENT VARIABLE",
                ]
                self.message(" ".join(msg))
                # TODO: Turn the sys.exit into exceptions.
                sys.exit(1)

        if (
            not filename
            or not os.path.isfile(filename)
            or not os.access(filename, os.R_OK)
        ):
            self.message(f"UNABLE TO READ: {filename}")
            # TODO: Turn the sys.exit into exceptions.
            sys.exit(2)

        self.__data: dict = dict()

        self._joined = None

        self.parsefile(
            filename,
            cache=keywords.get("cache", True),
            check=keywords.get("check", True),
        )

    def parsefile(self, filename: str, **keywords) -> None:
        """
        Method to parse or restore the database from cache.
        Called by :meth:`amespahdbpythonsuite.amespahdb.__init__` method.

        Parameters:
            filname : str
                String of the PAHdb filename and path.
            keywords
                Arbitrary keyword arguments.

        """

        # Create MD5 hash function pickle.
        md5 = (
            tempfile.gettempdir()
            + "/"
            + hashlib.md5(open(filename, "r").read().encode()).hexdigest()
            + ".pkl"
        )

        # Check if database is dumped in cache and restore it.
        if (
            keywords.get("cache", True)
            and os.path.isfile(md5)
            and os.access(md5, os.R_OK)
        ):
            self.message("RESTORING DATABASE FROM CACHE")

            # Start timer.
            tstart = time.perf_counter()

            # Open and read the dumped database.
            with open(md5, "rb") as f:
                self.__data = pickle.load(f)

            # Store the dumped database filename.
            self.__data["filename"] = md5

            # Stop timer and calculate elapsed time.
            elapsed = timedelta(seconds=(time.perf_counter() - tstart))

            info = [
                f"FILENAME                    : {md5}",
                f"ORIGNINAL FILENAME          : {filename}",
                f"PARSE TIME                  : {elapsed}",
                f'DATABASE                    : {self.__data["database"]}',
                f'VERSION (DATE)              : {self.__data["version"]} ({self.__data["date"]})',
                f'COMMENT                     : {self.__data["comment"]}',
            ]
            self.message(info, space=0)

        else:
            self.message("PARSING DATABASE: THIS MAY TAKE A FEW MINUTES")

            # Start timer.
            tstart = time.perf_counter()

            # Call XMLparser module to parse the database.
            parser = XMLparser(filename=filename, validate=keywords.get("check", True))

            # Store the database into self.__self.data.
            self.__data = parser.to_pahdb_dict()

            # Dump the database into pickle in the cache directory.
            if keywords.get("cache", True):
                with open(md5, "wb") as f:
                    pickle.dump(self.__data, f, pickle.HIGHEST_PROTOCOL)

            # Store the dumped database filename.
            self.__data["filename"] = filename

            # Stop timer and calculate elapsed time.
            elapsed = timedelta(seconds=(time.perf_counter() - tstart))

            info = [
                f'FILENAME                    : {self.__data["filename"]}',
                f"PARSE TIME                  : {elapsed}",
                f'DATABASE                    : {self.__data["database"]}',
                f'VERSION (DATE)              : {self.__data["version"]} ({self.__data["date"]})',
                f'COMMENT                     : {self.__data["comment"]}',
            ]

            self.message(info, space=0)

    def __repr__(self) -> str:
        """
        Class representation.

        """
        return f"{self.__class__.__name__}(" f"filename={self.__data['filename']})"

    def __str__(self) -> str:
        """
        A description of the instance.
        """

        return "AmesPAHdbPythonSuite AmesPAHdb instance."

    def __getkeybyuids(self, key: str, uids: list[int]) -> dict:
        """
        Get a dictionary of PAHdb properties
        retrieved by keyword for provided UIDs.

        Parameters:
            key : str
                Database keyword.
            uids : list of integers
                List of UIDs.

        Returns:
            Dictionary of retrieved properties with UIDs as keys.

        """

        if key == "species":
            return copy.deepcopy(
                dict(
                    (uid, self.__data["species"][uid])
                    for uid in uids
                    if uid in self.__data["species"].keys()
                )
            )
        else:
            return copy.deepcopy(
                dict(
                    (uid, self.__data["species"][uid][key])
                    for uid in uids
                    if uid in self.__data["species"].keys()
                )
            )

    def gettransitionsbyuid(self, uids: Union[list[int], int]) -> Transitions:
        """
        Retrieve and return transitions instance based on UIDs input.
        UIDs should be a list, e.g. the output of search method.
        Calls the :class:`amespahdbpythonsuite.transitions.Transitions` class.

        Parameters:
            uids : list of integers
                List of UIDs.

        Returns:
            transitions instance

        """

        uids_list = list()
        if isinstance(uids, int):
            uids_list.append(uids)
        else:
            uids_list = uids

        d = self.__getkeybyuids("transitions", uids_list)

        from amespahdbpythonsuite import transitions

        return transitions.Transitions(
            database=self.__data["database"],
            version=self.__data["version"],
            data=d,
            pahdb=self.__data,
            uids=list(d.keys()),
            model={"type": "zerokelvin_m", "temperature": 0.0, "description": ""},
            units={
                "abscissa": {
                    "unit": u.cm**-1,
                    "label": "frequency",
                },
                "ordinate": {"unit": u.km / u.mol, "label": "integrated cross-section"},
            },
        )

    def getlaboratorybyuid(self, uids: Union[list[int], int]) -> Optional[Laboratory]:
        """
        Retrieve and return laboratory database instance based on UIDs input.
        UIDs should be a list, e.g. the output of search method.
        Calls the :class:`amespahdbpythonsuite.laboratory.Laboratory` class.

        Parameters:
            uids : list of integers
                List of UIDs.

        Returns:
            laboratory database instance

        """

        # Check if the experimental database is loaded.
        if self.__data["database"] != "experimental":
            self.message("EXPERIMENTAL DATABASE REQUIRED")
            return None

        uids_list = list()
        if isinstance(uids, int):
            uids_list.append(uids)
        else:
            uids_list = uids

        d = self.__getkeybyuids("laboratory", uids_list)

        from amespahdbpythonsuite import laboratory

        return laboratory.Laboratory(
            database=self.__data["database"],
            version=self.__data["version"],
            data=d,
            pahdb=self.__data,
            uids=list(d.keys()),
            model={"type": "laboratory_m", "temperature": 0.0, "description": ""},
            units={
                "abscissa": {"unit": u.cm**-1, "label": "frequency"},
                "ordinate": {
                    "unit": u.def_unit(
                        "absorbance",
                        format={"latex": r"-\log(I/I_{0})"},
                        doc="Absorbance",
                    ),
                    "label": "absorbance",
                },
            },
        )

    def getspeciesbyuid(self, uids: Union[list[int], int]) -> Species:
        """
        Retrieve and return species instance based on UIDs input.
        UIDs should be a list, e.g. the output of search method.
        Calls the :class:`amespahdbpythonsuite.species.Species` class.

        Parameters:
            uids : list of integers
                List of UIDs.

        Returns:
            species instance

        """

        uids_list = list()
        if isinstance(uids, int):
            uids_list.append(uids)
        else:
            uids_list = uids

        d = self.__getkeybyuids("species", uids_list)

        from amespahdbpythonsuite import species

        return species.Species(
            database=self.__data["database"],
            version=self.__data["version"],
            data=d,
            pahdb=self.__data,
            uids=list(d.keys()),
            references=self.__getkeybyuids("references", uids_list),
            comments=self.__getkeybyuids("comments", uids_list),
        )

    def getgeometrybyuid(self, uids: Union[list[int], int]) -> Geometry:
        """
        Retrieve and return geometry instance based on UIDs input.
        UIDs should be a list, e.g. the output of search method.
        Calls the :class:`amespahdbpythonsuite.geometry.Geometry` class
        and :meth:`amespahdbpythonsuite.amespahdb.__getkeybyuids` method.

        Parameters:
            uids : list of integers
                List of UIDs.

        Returns:
            geometry instance

        """

        uids_list = list()
        if isinstance(uids, int):
            uids_list.append(uids)
        else:
            uids_list = uids

        d = self.__getkeybyuids("geometry", uids_list)

        from amespahdbpythonsuite import geometry

        return geometry.Geometry(
            database=self.__data["database"],
            version=self.__data["version"],
            data=d,
            pahdb=self.__data,
            uids=list(d.keys()),
        )

    def search(self, query: str) -> Optional[list]:
        """
        Search the database based on query input.

        Parameters:
            query : str
                String containing search query.

        Returns:
            List of UIDs

        Example:
            ``search('magnesium=0 oxygen=0 iron=0
            silicium=0 chx=0 ch2=0 c>20 h>0')``

        """

        if not query:
            return None

        words = list()

        n = len(query)

        i = 0
        while True:
            if i == n:
                break

            while query[i] == " ":
                i += 1

            token = query[i]
            if token in ["=", "<", ">", "(", ")"]:
                i += 1
                if i < n and query[i] == "=":
                    token += query[i]
                    i += 1
            elif token == "&":
                i += 1
                if i < n and query[i] == "&":
                    token += query[i]
                    i += 1
            elif token == "|":
                i += 1
                if i < n and query[i] == "|":
                    token += query[i]
                    i += 1
            elif token == "!":
                i += 1
                if i < n and query[i] == "=":
                    token += query[i]
                    i += 1
            else:
                i += 1
                while i < n and query[i] not in [
                    " ",
                    "=",
                    "<",
                    ">",
                    "&",
                    "|",
                    "(",
                    ")",
                    "!",
                ]:
                    token += query[i]
                    i += 1

            words.append(token)

        tokens = list()

        for word in words:
            tokens.append(self._tokenize(word))

        code = self._parsetokens(tokens)

        if not code:
            return None

        found = eval(
            f"[item[0] for item in _AmesPAHdb__data['species'].items() if ({code})]",
            self.__dict__ | {"np": np},
        )

        return found

    def _tokenize(self, word: str) -> dict:
        """
        A method called by :sec:`amespahdbpythonsuite.amespahdb.search`
        that creates a dictionary with keys based on the input word/category.

        Parameters:
            word : str

        """

        word = word.lower()

        token = {"type": "", "translation": "", "valid": False}

        charge = {
            "anion": 'item[1]["charge"] < 0',
            "cation": 'item[1]["charge"] > 0',
            "neutral": 'item[1]["charge"] == 0',
            "positive": 'item[1]["charge"] > 0',
            "negative": 'item[1]["charge"] < 0',
            "-": 'item[1]["charge"] == -1',
            "+": 'item[1]["charge"] == 1',
            "++": 'item[1]["charge"] == 2',
            "+++": 'item[1]["charge"] == 3',
            "---": 'item[1]["charge"] == -3',
        }

        identities = {
            "uid": "item[0]",
            "identifier": "item[0]",
            "hydrogen": 'item[1]["n_h"]',
            "carbon": 'item[1]["n_c"]',
            "nitrogen": 'item[1]["n_n"]',
            "oxygen": 'item[1]["n_o"]',
            "magnesium": 'item[1]["n_mg"]',
            "silicium": 'item[1]["n_si"]',
            "iron": 'item[1]["n_fe"]',
            "h": 'item[1]["n_h"]',
            "c": 'item[1]["n_c"]',
            "n": 'item[1]["n_n"]',
            "o": 'item[1]["n_o"]',
            "mg": 'item[1]["n_mg"]',
            "si": 'item[1]["n_si"]',
            "fe": 'item[1]["n_fe"]',
            "ch": 'item[1]["n_ch"]',
            "ch2": 'item[1]["n_ch2"]',
            "ch3": 'item[1]["n_ch3"]',
            "chx": 'item[1]["n_chx"]',
            "solo": 'item[1]["n_solo"]',
            "duo": 'item[1]["n_duo"]',
            "trio": 'item[1]["n_trio"]',
            "quartet": 'item[1]["n_quartet"]',
            "quintet": 'item[1]["n_quintet"]',
            "charge": 'item[1]["charge"]',
            "symmetry": 'item[1]["symmetry"]',
            "weight": 'item[1]["weight"]',
            "scale": 'item[1]["scale"]',
            "energy": 'item[1]["total_e"]',
            "zeropoint": 'item[1]["vib_e"]',
            "experiment": 'item[1]["exp"]',
        }

        strings: list[str] = []
        if self.__data["database"] == "clusters/theoretical":
            identities.update(
                {
                    "monomers": 'item[1]["monomers"]',
                    "type": 'item[1]["type"]',
                    "conformation": 'item[1]["conformation"]',
                }
            )
            strings.extend(
                set([s["monomers"] for s in self.__data["species"].values()])
            )
            strings.extend(set([s["type"] for s in self.__data["species"].values()]))
            strings.extend(
                set([s["conformation"] for s in self.__data["species"].values()])
            )

        composed = {
            "wavenumber": 'np.array([t["frequency"] {operator} {operand} for t in item[1]["transitions"]])',
            "absorbance": 'np.array([t["intensity"] {operator} {operand} for t in item[1]["transitions"]])',
            "frequency": 'np.array([t["frequency"] {operator} {operand} for t in item[1]["transitions"]])',
            "intensity": 'np.array([t["intensity"] {operator} {operand} for t in item[1]["transitions"]])',
        }

        logical = {"and": "and", "or": "or", "|": "or", "&": "and"}

        comparison = {
            "<": "<",
            "lt": "<",
            ">": ">",
            "gt": ">",
            "=": "==",
            "eq": "==",
            "<=": "<=",
            "le": "<=",
            ">=": ">=",
            "ge": ">=",
            "with": "and",
            "ne": "!=",
            "!=": "!=",
        }

        transfer = {"(": "(", ")": ")"}

        if word.strip("+-").isnumeric():
            token["type"] = "NUMERIC"
            token["translation"] = word
        elif word in charge:
            token["type"] = "CHARGE"
            token["translation"] = charge[word]
        elif word in identities:
            token["type"] = "IDENTITY"
            token["translation"] = identities[word]
        elif word in composed:
            token["type"] = "COMPOSED"
            token["translation"] = composed[word]
        elif word in logical:
            token["type"] = "LOGICAL"
            token["translation"] = logical[word]
        elif word in comparison:
            token["type"] = "COMPARISON"
            token["translation"] = comparison[word]
        elif word in transfer:
            token["type"] = "TRANSFER"
            token["translation"] = transfer[word]
        elif word in strings:
            token["type"] = "STRING"
            token["translation"] = f"'{word}'"
        elif re.search(
            "(mg+|si+|fe+|[chno]+)([0-9]*)(mg+|si+|fe+|[chno]+)([0-9]*)(mg+|si+|fe+|[chno]*)([0-9]*)",
            word,
        ):
            token["type"] = "FORMULA"
            token["translation"] = f"item[1]['formula'] == '{word.upper()}'"
        else:
            # TODO: add search by compound name
            token["type"] = "IGNORE"
            token["translation"] = word
            return token

        token["valid"] = True

        return token

    def _parsetokens(self, tokens: list) -> str:
        """
        Parse the dictionary of tokens created by
        :sec:`amespahdbpythonsuite.amespahdb.tokenize`
        and return string of expressions.

        Parameters:
            tokens : list
                List of dictionaries.

        Returns:
            parsed : str
                String of expressions based on tokens.

        """

        ntokens = len(tokens)

        prev = -1

        current = 0

        next = 1 if ntokens > 1 else -1

        parsed = ""

        sub = ""

        while current != -1:
            if tokens[current]["type"] == "FORMULA":
                if prev > -1:
                    if not (
                        tokens[prev]["type"] != "LOGICAL" and tokens[prev]["valid"]
                    ):
                        parsed += " or "
                parsed += tokens[current]["translation"]
            elif tokens[current]["type"] == "IDENTITY":
                if prev > -1:
                    if not (
                        tokens[prev]["type"] == "LOGICAL"
                        or tokens[prev]["type"] == "TRANSFER"
                        and tokens[prev]["valid"]
                    ):
                        parsed += " and "
                if next > -1:
                    if tokens[next]["type"] == "COMPARISON":
                        parsed += " " + tokens[current]["translation"]
                    else:
                        parsed += " " + tokens[current]["translation"] + " > 0"
            elif tokens[current]["type"] == "COMPOSED":
                if prev > -1:
                    if not (
                        tokens[prev]["type"] == "LOGICAL"
                        or tokens[prev]["type"] == "TRANSFER"
                        and tokens[prev]["valid"]
                    ):
                        if sub:
                            sub += " & "
                        else:
                            parsed += " and "

                if next > -1:
                    if tokens[next]["type"] != "COMPARISON":
                        print("EXPECTING OPERATOR")
                        return ""

                    partial = tokens[current]["translation"].replace(
                        "{operator}", tokens[next]["translation"]
                    )

                    prev = current

                    current = next

                    if next:
                        if next == ntokens - 1:
                            next = -1
                        else:
                            next += 1

                    if next > -1 and tokens[next]["type"] != "NUMERIC":
                        print("EXPECTING OPERAND")
                        return ""

                    sub += " " + partial.replace(
                        "{operand}", tokens[next]["translation"]
                    )

                    prev = current

                    current = next

                    if next:
                        if next == ntokens - 1:
                            next = -1
                        else:
                            next += 1
            elif tokens[current]["type"] == "NUMERIC":
                if prev > -1:
                    if tokens[prev]["type"] == "COMPARISON" and tokens[prev]["valid"]:
                        parsed += " " + tokens[current]["translation"]
                    else:
                        tokens[current]["valid"] = False
            elif tokens[current]["type"] == "LOGICAL":
                if prev > -1:
                    if (
                        tokens[prev]["type"] == "IDENTITY"
                        or tokens[prev]["type"] == "NUMERIC"
                        or tokens[prev]["type"] == "FORMULA"
                        or tokens[prev]["type"] == "CHARGE"
                        and tokens[prev]["valid"]
                    ):
                        if next > -1:
                            if tokens[next]["type"] == "TRANSFER":
                                parsed += tokens[current]["translation"]
                            elif tokens[next]["type"] == "COMPOSED" and sub:
                                sub += {"and": " & ", "or": " | "}[
                                    tokens[current]["translation"]
                                ]
                            elif (
                                tokens[next]["type"] == "IDENTITY"
                                or tokens[next]["type"] == "NUMERIC"
                                or tokens[next]["type"] == "FORMULA"
                                or tokens[next]["type"] == "CHARGE"
                            ):
                                parsed += " " + tokens[current]["translation"]
                            else:
                                tokens[current]["valid"] = False

            elif tokens[current]["type"] == "COMPARISON":
                if prev > -1:
                    if tokens[prev]["type"] == "IDENTITY" and tokens[prev]["valid"]:
                        if next is not None:
                            if tokens[next]["type"] in ["NUMERIC", "STRING"]:
                                parsed += " " + tokens[current]["translation"]
                            else:
                                tokens[current]["valid"] = False
            elif tokens[current]["type"] == "CHARGE":
                if prev > -1:
                    if not (
                        tokens[prev]["type"] == "LOGICAL" and tokens[prev]["valid"]
                    ):
                        parsed += " and "
                parsed += " " + tokens[current]["translation"]
            elif tokens[current]["type"] == "TRANSFER":
                parsed += tokens[current]["translation"]
            elif tokens[current]["type"] == "NAME":
                if prev > -1:
                    if not (
                        tokens[prev]["type"] == "LOGICAL" and tokens[prev]["valid"]
                    ):
                        parsed += " and "
                # TODO implement name
                # parsed += f"item[1]['comments'] == {tokens[current]['translation']}"
            elif tokens[current]["type"] == "STRING":
                if prev > -1:
                    if tokens[prev]["type"] in ["COMPARISON", "NUMERIC"] and tokens[prev]["valid"]:
                        parsed += tokens[current]["translation"]
                    else:
                        print("EXPECTING COMPARISON")
                        tokens[current]["valid"] = False
                else:
                    print("EXPECTING COMPARISON")
                    tokens[current]["valid"] = False

            elif tokens[current]["type"] == "IGNORE":
                print(f"'{tokens[current]['translation']}' NOT UNDERSTOOD")

                return ""

            prev = current

            current = next

            if next:
                if next == ntokens - 1:
                    next = -1
                else:
                    next += 1

        if sub:
            parsed += f" np.any({sub})"

        return parsed

    def getversion(self) -> str:
        """
        Method to retrieve the PAHdb version.

        Returns:
            String of PAHdb version.

        """
        return self.__data["version"]

    def checkversion(self, version: str) -> bool:
        """
        Method to check against a PAHdb version.

        Returns:
            Boolean whether a provided version matched the PAHdb version.

        """
        return version == self.__data["version"]

    def gettype(self) -> str:
        """
        Method to retrieve the PAHdb type.

        Returns:
            String of PAHdb type.

        """
        return self.__data["database"]

    def getdatabaseref(self) -> dict:
        """
        Method to retrieve the database.

        Returns:
            Dictionary containing the parsed database.

        """
        return self.__data

    def extenddatabase(self, **keywords) -> None:
        """
        Method to extend the database.

        Returns:
            None

        """
        if not self.__data:
            self.message("NO VALID DATABASE")

        if keywords.get("species", False):
            self.__data["species"].update(keywords["species"])

    @staticmethod
    def message(text, space: int = 55) -> None:
        """
        A method to print terminal message.

        Parameters:
            text : string or list of strings.
                Text to be displayed.
            space : integer
                Number to indent the text.

        """
        line = (space + 2) * "="
        print(line)
        if isinstance(text, list):
            for t in text:
                print(t.center(space))
        else:
            print(text.center(space))
        print(line)
        print()
