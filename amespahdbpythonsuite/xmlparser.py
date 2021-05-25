#!/usr/bin/env python3
"""xmlparser.py

Parse a NASA Ames PAH IR Spectroscopic Database XML library file, with
or without schema checking.

"""

import array
import base64
import urllib.request

from lxml import etree


class XMLparser:
    """Parse a NASA Ames PAH IR Spectroscopic library XML-file.

    Optional behaviour includes validating against a schema.

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

    def __init__(self, filename=None, validate=False):
        """Inits XMLparser with schema checking off, no given filename."""
        self.filename = filename
        self.validate = validate
        self.library = {}

    def __repr__(self):
        """Class representation."""
        return (f'{self.__class__.__name__}('
                f'filename={self.filename!r}, ')

    def __str__(self):
        """A description of the object."""
        return (f'XML parser from the AmesPAHdbPythonSuite. '
                f'XML file: \n{self.filename!r}.')

    def verify_schema(self):
        """Validate against linked schema.

        Note:
            Requires that self.filename is set.

            It sets the internal attributes:
            _tree (etree.ElementTree): parsed XML file.
            _root: root element of ElementTree tree.

        Returns:
            True if successful, otherwisse False.

        """

        self._tree = etree.parse(self.filename)
        self._root = self._tree.getroot()

        schema = self._root.get(
            '{http://www.w3.org/2001/XMLSchema-instance}' + 'schemaLocation'
        )

        if schema:
            _, uri = schema.split(' ', 1)
            doc = etree.parse(urllib.request.urlopen(uri))
            xmlschema = etree.XMLSchema(doc)
            try:
                xmlschema.assertValid(self._tree)
            except Exception as e:
                raise e
            else:
                return True

    def to_pahdb_dict(self, validate=False):
        """Parses the XML, with or without validation.

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

            self._context = \
                etree.iterwalk(self._tree, events=("start", "end"))
        else:
            self._context = \
                etree.iterparse(self.filename, events=("start", "end"))

        self.library = self._tree_to_pahdb_dict()

        return self.library

    def _tree_to_pahdb_dict(self):
        """Convert the element tree to a a pahdb_dict.

        Returns: library: Dictionary, with the UIDs as keys,
            containing the transitions, geometry data, as well as UID
            metadata, references, comments, and laboratory.

        """

        while True:
            action, elem = next(self._context)
            tag = etree.QName(elem).localname

            if action == 'start':
                if tag == 'species':
                    self.library['species'] = \
                        self._species_handler(self._context)
                elif tag == 'pahdatabase':
                    self.library.update(elem.attrib)
            elif action == 'end':
                if tag == 'comment':
                    self.library['comment'] = elem.text
                elif tag == 'pahdatabase':
                    break
                elem.clear()
        return self.library

    def _species_handler(self, context):
        """Parse a PAHdb XML <species> tag."""
        species = {}

        while True:
            action, elem = next(context)
            tag = etree.QName(elem).localname

            if action == 'start' and tag == 'specie':
                uid = int(elem.attrib['uid'])
                species[uid] = self._specie_handler(context)
            elif action != 'end':
                continue

            if tag == 'species':
                break

            elem.clear()

        return species

    def _specie_handler(self, context):
        """Parse a PAHdb XML <specie> tag."""

        def specie_geometry_handler(context):
            """<specie> tag: Parse its child <geometry> tag."""
            geometry = []

            while True:
                action, elem = next(context)
                tag = etree.QName(elem).localname

                if tag == 'atom' and action == 'start':
                    atom_dict = {}

                    while True:
                        action, elem = next(context)
                        tag = etree.QName(elem).localname

                        if action != 'end':
                            continue

                        if tag == 'atom':
                            geometry.append(atom_dict)
                            break

                        atom_dict[tag] = float(elem.text)

                        elem.clear()
                        
                elif action != 'end':
                    continue

                if tag == 'geometry':
                    break

                elem.clear()

            return geometry

        def specie_transitions_handler(context):
            """<specie> tag: Parse its child <transitions> tag."""
            transitions = []

            while True:
                action, elem = next(context)
                tag = etree.QName(elem).localname

                if tag == 'mode':
                    mode_dict = {}

                    while True:
                        action, elem = next(context)
                        tag = etree.QName(elem).localname

                        if action != 'end':
                            continue

                        if tag == 'mode':
                            transitions.append(mode_dict)
                            break

                        if elem.attrib:
                            mode_dict.update(elem.attrib)

                        try:
                            value = float(elem.text)
                        except ValueError:
                            value = elem.text
                        mode_dict[tag] = value

                        elem.clear()
                elif action != 'end':
                    continue

                if tag == 'transitions':
                    break

                elem.clear()

            return transitions

        def specie_laboratory_handler(context):
            """<specie> tag: Parse its child <laboratory> tag."""
            laboratory = {}

            while True:
                action, elem = next(context)
                tag = etree.QName(elem).localname

                if action == 'end':
                    if tag == 'frequency' or tag == 'intensity':
                        bin = base64.b64decode(elem.text)
                        laboratory[tag] = array.array('f', bin)
                    elif tag == 'laboratory':
                        break

                    elem.clear()

            return laboratory

        specie_dict = {}

        while True:
            action, elem = next(context)
            tag = etree.QName(elem).localname

            if action == 'start':
                if tag == 'comments':
                    comments = []

                    while True:
                        action, elem = next(context)
                        tag = etree.QName(elem).localname

                        if action == 'end':
                            if tag == 'comment':
                                comments.append(elem.text)
                            elif tag == 'comments':
                                break

                            elem.clear()

                    specie_dict['comments'] = tuple(comments)
                elif tag == 'references':
                    references = []

                    while True:
                        action, elem = next(context)
                        tag = etree.QName(elem).localname

                        if action == 'end':
                            if tag == 'reference':
                                references.append(elem.text)
                            elif tag == 'references':
                                break

                            elem.clear()

                    specie_dict['references'] = tuple(references)
                elif tag == 'geometry':
                    specie_dict['geometry'] = specie_geometry_handler(context)
                elif tag == 'transitions':
                    specie_dict['transitions'] = \
                        specie_transitions_handler(context)
                elif tag == 'laboratory':
                    specie_dict['laboratory'] = \
                        specie_laboratory_handler(context)
            elif action == 'end':
                if tag == 'specie':
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
