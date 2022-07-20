#!/usr/bin/env python3
"""xmlparser.py

Parse a NASA Ames PAH IR Spectroscopic Database XML library file, with
or without schema checking.

"""

from typing import Union
import array
import base64
import urllib.request
from urllib.error import HTTPError, URLError

from lxml import etree  # type: ignore


class XMLparser:
    """
    Parse a NASA Ames PAH IR Spectroscopic library XML-file.

    Optional behavior includes validating against a schema.

    Attributes:
        filename (str): XML filename.
        validate (bool): Whether to validate the XML schema or not.
        library (dict): Dictionary containing the parsed data.

    Examples:
        parser = XMLparser(filename="xml_file.xml")
        parser.verify_schema()
        dict = parser.to_pahdb_dict()

        parser = XMLparser(filename="xml_file.xml", validate=True)
        dict = parser.to_pahdb_dict()

        parser = XMLparser()
        parser.filename = "xml_file.xml"
        parser.validate = True
        dict = parser.to_pahdb_dict()

    """

    def __init__(self, filename: str = None, validate: bool = False) -> None:
        """
        Inits XMLparser with schema checking off, no given filename.

        """
        self.filename = filename
        self.validate = validate
        self.library: dict = dict()

    def __repr__(self) -> str:
        """
        Class representation.

        """
        return f"{self.__class__.__name__}(" f"{self.filename=})"

    def __str__(self) -> str:
        """
        A description of the instance.

        """
        return f"AmesPAHdbPythonSuite XML parser.\n" f"XML filename: {self.filename=}."

    def verify_schema(self) -> bool:
        """
        Validate against linked schema.

        Note:
            Requires that self.filename is set.

            It sets the internal attributes:
            _tree (etree.ElementTree): parsed XML file.
            _root: root element of ElementTree tree.

        Returns:
            True if successful, otherwise False.

        """

        self._tree = etree.parse(self.filename)
        self._root = self._tree.getroot()

        schema = self._root.get(
            "{http://www.w3.org/2001/XMLSchema-instance}" + "schemaLocation"
        )

        if schema:
            _, uri = schema.split(" ", 1)
            try:
                response = urllib.request.urlopen(uri, timeout=3.0)
            except (HTTPError, URLError):
                # TODO For now, fallback to True if we can't get a schema, use False instead?
                return True
            doc = etree.parse(response)
            xmlschema = etree.XMLSchema(doc)
            try:
                xmlschema.assertValid(self._tree)
            except Exception as e:
                raise e
            else:
                return True

        return True

    def to_pahdb_dict(self, validate: bool = False) -> dict:
        """
        Parses the XML, with or without validation.

        Args:
            validate (bool). Defaults to self.valdiate value, but can be
            overridden.

        Note:
            Sets the attribute self.library when successful.

        Returns: library (dict): Dictionary, with the UIDs as keys,
            containing the transitions, geometry data, as well as UID
            metadata, references, comments, and laboratory.

        """

        if self.validate or validate:
            self.verify_schema()
            self._context = etree.iterwalk(self._tree, events=("start", "end"))
        else:
            self._context = etree.iterparse(self.filename, events=("start", "end"))

        self.library = self._tree_to_pahdb_dict()

        return self.library

    def _tree_to_pahdb_dict(self) -> dict:
        """
        Convert the element tree to a a pahdb_dict.

        Returns: library: Dictionary, with the UIDs as keys,
            containing the transitions, geometry data, as well as UID
            metadata, references, comments, and laboratory.

        """

        while True:
            action, elem = next(self._context)
            tag = etree.QName(elem).localname

            if action == "start":
                if tag == "species":
                    self.library["species"] = self._species_handler(self._context)
                elif tag == "pahdatabase":
                    self.library.update(elem.attrib)
            elif action == "end":
                if tag == "comment":
                    self.library["comment"] = elem.text
                elif tag == "pahdatabase":
                    break
                elem.clear()
        return self.library

    def _species_handler(self, context: Union[etree.iterwalk, etree.iterparse]) -> dict:
        """
        Parse a PAHdb XML <species> tag.

        """
        species = dict()

        while True:
            action, elem = next(context)
            tag = etree.QName(elem).localname

            if action == "start" and tag == "specie":
                uid = int(elem.attrib["uid"])
                species[uid] = self._specie_handler(context)
            elif action != "end":
                continue

            if tag == "species":
                break

            elem.clear()

        return species

    def _specie_handler(self, context: Union[etree.iterwalk, etree.iterparse]) -> dict:
        """
        Parse a PAHdb XML <specie> tag.

        """

        def specie_geometry_handler(
            context: Union[etree.iterwalk, etree.iterparse]
        ) -> list:
            """<specie> tag: Parse its child <geometry> tag."""
            geometry = list()

            while True:
                action, elem = next(context)
                tag = etree.QName(elem).localname

                if tag == "atom" and action == "start":
                    atom_dict: dict = dict()

                    while True:
                        action, elem = next(context)
                        tag = etree.QName(elem).localname

                        if action != "end":
                            continue

                        if tag == "atom":
                            geometry.append(atom_dict)
                            break

                        atom_dict[tag] = float(elem.text)

                        elem.clear()

                elif action != "end":
                    continue

                if tag == "geometry":
                    break

                elem.clear()

            return geometry

        def specie_transitions_handler(
            context: Union[etree.iterwalk, etree.iterparse]
        ) -> list:
            """
            <specie> tag: Parse its child <transitions> tag.

            """
            transitions = list()

            while True:
                action, elem = next(context)
                tag = etree.QName(elem).localname

                if tag == "mode":
                    mode_dict: dict = dict()

                    while True:
                        action, elem = next(context)
                        tag = etree.QName(elem).localname

                        if action != "end":
                            continue

                        if tag == "mode":
                            transitions.append(mode_dict)
                            break

                        if elem.attrib:
                            for attr, text in elem.attrib.items():
                                try:
                                    value = float(text)
                                except ValueError:
                                    value = text
                                mode_dict[attr] = value

                        try:
                            value = float(elem.text)
                        except ValueError:
                            value = elem.text
                        mode_dict[tag] = value

                        elem.clear()
                elif action != "end":
                    continue

                if tag == "transitions":
                    break

                elem.clear()

            return transitions

        def specie_laboratory_handler(
            context: Union[etree.iterwalk, etree.iterparse]
        ) -> dict:
            """
            <specie> tag: Parse its child <laboratory> tag.

            """
            laboratory = dict()

            while True:
                action, elem = next(context)
                tag = etree.QName(elem).localname

                if action == "end":
                    if tag == "frequency" or tag == "intensity":
                        bin = base64.b64decode(elem.text)
                        laboratory[tag] = array.array("f", bin)
                    elif tag == "laboratory":
                        break

                    elem.clear()

            return laboratory

        specie_dict: dict = {
            "comments": list(),
            "references": list(),
            "geometry": list(),
            "transitions": list(),
            "laboratory": dict(),
        }

        while True:
            action, elem = next(context)
            tag = etree.QName(elem).localname

            if action == "start":
                if tag == "comments":
                    comments = list()

                    while True:
                        action, elem = next(context)
                        tag = etree.QName(elem).localname

                        if action == "end":
                            if tag == "comment":
                                comments.append(elem.text)
                            elif tag == "comments":
                                break

                            elem.clear()

                    specie_dict["comments"] = comments
                elif tag == "references":
                    references = list()

                    while True:
                        action, elem = next(context)
                        tag = etree.QName(elem).localname

                        if action == "end":
                            if tag == "reference":
                                references.append(elem.text)
                            elif tag == "references":
                                break

                            elem.clear()

                    specie_dict["references"] = references
                elif tag == "geometry":
                    specie_dict["geometry"] = specie_geometry_handler(context)
                elif tag == "transitions":
                    specie_dict["transitions"] = specie_transitions_handler(context)
                elif tag == "laboratory":
                    specie_dict["laboratory"] = specie_laboratory_handler(context)
            elif action == "end":
                if tag == "specie":
                    break

                try:
                    value = float(elem.text)
                    if value % 1 == 0:
                        value = int(value)
                except ValueError:
                    value = elem.text
                specie_dict[tag] = value

                elem.clear()

        return specie_dict
