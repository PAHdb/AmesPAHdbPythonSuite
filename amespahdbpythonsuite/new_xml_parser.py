#!/usr/bin/env python3
"""
parsers.py

Parse a PAHdb XML file, with or without schema checking.
"""

import array
import base64
import os
import urllib.request

import numpy as np
from lxml import etree as ElementTree

from amespahdbpythonsuite.utils.xml_data_types import \
    metadata_keys, geometry_keys, transition_keys, laboratory_keys

from amespahdbpythonsuite.handlers import pahdatabase_handler


from ipdb import set_trace as st


class parsey:

    def __init__(self, filename):
        self.filename = filename

        self._tree = ElementTree.parse(self.filename)
        self._root = self._tree.getroot()
        self._schema = self._root.get(
            '{http://www.w3.org/2001/XMLSchema-instance}' + 'schemaLocation'
        )

        # self.parse()

    def parse(self):
        print("hi")
        st()













class XMLparser:
    """Parse a PAHdb XML file.

    Optional behaviour includes validating the XML schema.

    Attributes:
        filename (str): XML file path.
        validate (bool): Whether to validate the XML schema or not.
        schema_is_valid (bool): Outcome of the schema checking.
        database (dict): Dictionary containing the PAHdb-parsed data.
        info (dict): Dictionary containing the PAHdb-parsed metadata.

    Note:
        `validated` (bool) indicates if the schema is valid (True) or
            not (False). Set to None until validation is performed.

    Examples:
        parser = XMLparser(filename="xml_file.xml")
        parser.verify_schema()
        database, info = parser.to_pahdb_dict()

        parser = XMLparser(filename="xml_file.xml", validate=True)
        database, info = parser.to_pahdb_dict()

        parser = XMLparser()
        parser.filename = "xml_file.xml"
        parser.validate = True
        database, info = parser.to_pahdb_dict()

        parser = XMLparser.from_string("<XML code here>")
        parser.verify_schema()
        database, info = parser.to_pahdb_dict()

    """

    def __init__(self, filename=None, validate=False):
        """Inits XMLparser with schema checking off, no given filename."""
        self.filename = filename
        self.validate = validate
        self.schema_is_valid = None

        self.database = None
        self.info = None

    def __repr__(self):
        """Class representation."""
        return (f'{self.__class__.__name__}('
                f'filename={self.filename!r}, '
                f'validate={self.validate!r})')

    def __str__(self):
        """A description of the object."""
        return (f'An XML parser from the AmesPAHdbPythonSuite. '
                f'XML file: \n{self.filename!r}.')

    def verify_schema(self):
        """Validate XML schema.

        Note:
            Requires that self.filename is set.

            It sets the internal attributes:
            _tree (ElementTree.ElementTree): parsed XML file.
            _root: root element of ElementTree tree.

        Returns:
            True if successful. Also assigns self.schema_is_valid to True
            if successful, False if unsuccessful.
        """

        self._tree = ElementTree.parse(self.filename)
        self._root = self._tree.getroot()

        schema = self._root.get(
            '{http://www.w3.org/2001/XMLSchema-instance}' + 'schemaLocation'
        )

        if schema:
            _, uri = schema.split(' ', 1)
            doc = ElementTree.parse(urllib.request.urlopen(uri))
            xmlschema = ElementTree.XMLSchema(doc)
            try:
                xmlschema.assertValid(self._tree)
            except Exception as e:
                self.schema_is_valid = False
                print("Schema invalid!")
                raise e
            else:
                self.schema_is_valid = True
                print("Successfully validated schema.")
                return True

    def to_pahdb_dict(self, validate=False):
        """Perform the actual parsing, with or without validation.

        Args:
            validate (bool). Defaults to self.valdiate value, but can be
            overridden.

        Note:
            Sets the attributes self.database and self.info if successful.

        Returns:
            database (dict): Dictionary, with the UIDs as keys,
                containing the transitions, geometry data, as well as
                UID metadata, references, comments.
            info (dict): Dictionary containing general information about
                the XML file, including filename, type, date, full
                (bool), version, comment, and number of species.
        """
        events = ("start", "end")

        # First determine whether we've already validated the schema.
        if self.schema_is_valid:
            self._formatted_tree = \
                ElementTree.iterwalk(self._tree, events=events)

        # Else decide whether to validate here, or skip validation.
        else:
            if self.validate or validate:
                self.verify_schema()
                self._formatted_tree = \
                    ElementTree.iterwalk(self._tree, events=events)
            else:
                self._formatted_tree = \
                    ElementTree.iterparse(self.filename, events=events)

        self.database, self.info = \
            self._tree_to_pahdb_dict(self._formatted_tree)

        return self.database, self.info

    def _tree_to_pahdb_dict(self, formatted_tree):
        """Convert the tree to a dict of dicts.

        Args:
            formatted_tree: lxml tree of PAHdb file.

        Returns:
            database_dict: Dictionary, with the UIDs as keys,
                containing the transitions, geometry data, as well as
                UID metadata, references, comments.
            info_dict: Dictionary containing general information about
                the XML file, including filename, type, date, full
                (bool), version, comment, and number of species.
        """
        # Advance to root and identify namespace.
        event, root = next(formatted_tree)
        namespace = self.determine_xml_namespace(root)

        # Root dict attributes.
        root_dict = root.attrib

        # Form important dicts, to be returned.
        database_dict = {}
        info_dict = {}
        first_run = True

        # Iterate through the tree, clearing elements as we go.
        for _, (event, elem) in enumerate(formatted_tree):

            # Store the top XML comment.
            if first_run:
                if event == "end" and elem.tag == namespace + 'comment':
                    pahdb_comment = elem.text
                    first_run = False

            # Extract UID data on a per-tag ('specie') basis.
            if event == "end" and elem.tag == namespace + 'specie':
                # Extract the metadata.
                uid, metadata, comments, reference = \
                    self._extract_info(elem, namespace)

                # Extract the numeric data.
                args = elem, namespace
                geometry_dict = self._extract_numeric(*args, 'geometry')
                transitions_dict = self._extract_numeric(*args, 'transitions')
                laboratory_dict = self._extract_numeric(*args, 'laboratory')

                # Construct the UID-specific dictionary.
                database_dict[uid] = {
                    'metadata': metadata,
                    'geometry': geometry_dict,
                    'transitions': transitions_dict,
                    'comments': comments,
                    'reference': reference,
                    'laboratory': laboratory_dict,
                }

                elem.clear()

        # Path of XML file for saving in dict.
        full_path = os.path.abspath(self.filename)

        # Sanity check, no repeated UIDs.
        db_keys = [key for key in database_dict]
        assert len(db_keys) == len(np.unique(db_keys))

        # Construct dict with general XML information.
        info_dict = {
            'filename': full_path,
            'type': root_dict['database'],
            'date': root_dict['date'],
            'full': bool(root_dict['full']),
            'version': root_dict['version'],
            'comment': pahdb_comment,
            'nspecies': len(db_keys),
        }

        return database_dict, info_dict

    @staticmethod
    def determine_xml_namespace(root):
        """Returns the namespace for an ElementTree root object.

        Args:
            root (ElementTree): root object.

        Returns:
            namespace (str):
        """
        namespace = ''
        try:
            namespace_end_index = root.tag.index('}')
        except ValueError:
            pass
        else:
            namespace = root.tag[:namespace_end_index + 1]

        return namespace

    @staticmethod
    def _extract_info(specie, namespace):
        """Returns a dictionary of PAH metadata for a single 'specie'.

        Args:
            specie: XML <specie> element.
            namespace (str): XML namespace.

        Returns:
            uid (int): Unique ID specific to PAHdb.
            sub_dict (dict): Contains all UID-specific info.
            comment_dict (dict): General comments from XML file.
            ref_list (tuple): Tuple of general references from XML file.
        """

        def extract_metadata(specie, namespace):
            """Returns the specie metadata.

            Args:
                specie: XML <specie> element.
                namespace (str): XML namespace.

            Returns:
                sub_dict (dict): Contains all UID-specific info.
                comment_dict (dict): General comments from XML file.
                ref_list (tuple): Tuple of general references from XML file.
            """

            # Create UID-specific dictionary, empty ref list.
            sub_dict = {}
            ref_list = []

            # Identify its children and their tags.
            children = specie.getchildren()
            childrens_tags = [x.tag for x in children]

            # # Retrieve metadata from each child.
            # metadata_keys = get_metadata_keys()

            # Iterative over the metadata keys.
            for (meta_key, meta_type) in list(metadata_keys.items()):

                # Identify location of metadata for meta_key.
                try:
                    metadata_index = \
                        childrens_tags.index(namespace + meta_key)
                except ValueError:
                    continue
                else:
                    # Isolate metadata.
                    metadata = children[metadata_index]

                # Additional logic for references tag, if present.
                if meta_key == 'references':
                    ref_children = metadata.getchildren()
                    ref_list = tuple([item.text for item in ref_children])
                    continue

                # Additional logic for (maybe) multi-element comments tag.
                if meta_key == 'comments':

                    # Iterate over list of comments, parse individually.
                    comment_list = []
                    for _, entry in enumerate(metadata.getchildren()):
                        assert entry.tag == namespace + 'comment'

                        # key-value pair for comment text.
                        tmp_dict = {
                            'comment': entry.text,
                        }

                        # key-value pair for comment type.
                        tmp_dict.update(entry.attrib)

                        # store each comment (dict) in the list.
                        comment_list.append(tmp_dict)

                    continue

                # Extract metadata content.
                metadata_content = metadata.text

                # Cast to desired type (key dependent).
                if not isinstance(metadata_content, meta_type):
                    metadata_content = meta_type(metadata_content)

                # Insert metadata into UID-specific dictionary.
                sub_dict[meta_key] = metadata_content

            return sub_dict, comment_list, ref_list

        # Unique identifier.
        try:
            uid = int(specie.attrib['uid'])
        except KeyError:
            uid = 0

        sub_dict, comment_dict, ref_list = \
            extract_metadata(specie, namespace)

        return uid, sub_dict, comment_dict, ref_list

    @staticmethod
    def _extract_numeric(specie, namespace, label):
        """Returns a dictionary (or tuple thereof) for numeric data
        from an XML 'specie'.

        Args:
            specie: XML <specie> element.
            namespace (str): XML namespace.
            label (str): to identify the element as 'geometry', 'transition',
                or 'laboratory' data.

        Returns:
            Tuple of dictionaries. The number of dicts is dependent on
            what type of data is trying to be extracted: if
            label == 'laboratory', a single dict containing the
            frequency/intensity information is present within the
            tuple. If label is in ['geometry', 'transitions'], one dict
            per atom in the geometry or vibrational mode is present in
            the tuple (respectively).

        Note:
            If laboratory data, returns a single dictionary with keys
            `frequency`, `intensity`. If geometry data, returns a tuple of
            geometry dictionaries (one per atom), each with the keys `x`, `y`,
            `z`, `position`. If transition data, returns a tuple of transition
            dictionaries (one per transition/mode), each with keys `frequency`,
            `scale`, `intensity`, `symmetry`.
        """

        def parse_atom_to_dict(atom, namespace, geometry_keys):
            """Parse a single <atom> XML element to dictionary.

            Args:
                atom: XML element of type <atom> to be parsed.
                namespace (str): XML namespace.
                geometry_keys (dict): Dict of expected keys/dtypes.

            Returns:
                atom_dict (dict): Dict of atom properties/quantities.
            """
            atom_properties = atom.getchildren()
            atom_dict = {}

            for atom_property in atom_properties:
                # Parse children.
                key = atom_property.tag.split(namespace)[-1]
                value = atom_property.text

                if key == 'type':
                    value = int(value)
                else:
                    # Cast to desired type (key dependent).
                    desired_type = geometry_keys[key]
                    if not isinstance(value, desired_type):
                        value = desired_type(value)

                # Store in new dictionary.
                atom_dict[key] = value

            return atom_dict

        def parse_mode_to_dict(mode, namespace, transition_keys):
            """Parse a single <mode> XML element to dictionary.

            Args:
                mode: XML element of type <mode> to be parsed.
                namespace (str): XML namespace.
                transition_keys (dict): Dict of expected keys/dtypes.

            Returns:
                mode_dict (dict): Dict of mode properties/quantities.
            """
            mode_properties = mode.getchildren()
            mode_dict = {}

            for mode_property in mode_properties:
                key = mode_property.tag.split(namespace)[-1]
                value = mode_property.text  # of type str.

                # Cast to desired type (key dependent).
                desired_type = transition_keys[key]
                if not isinstance(value, desired_type):
                    value = desired_type(value)

                # Store in new dictionary.
                mode_dict[key] = value

                # Extract from text field, just for 'frequency'.
                if key == 'frequency':
                    attrib = mode_property.attrib  # of type dict.
                    if 'scale' in attrib:  # theoretical data XML file.
                        assert len(list(attrib)) == 1
                        assert list(attrib)[0] == 'scale'
                        scale_factor = attrib['scale']
                        mode_dict['scale'] = float(scale_factor)

            return mode_dict

        def parse_lab_to_dict(data, namespace, lab_keys):
            """Parse a single <laboratory> XML element to dictionary.

            Args:
                data: XML element of type <laboratory> to be parsed.
                namespace (str): XML namespace.
                lab_keys (dict): Dict of expected keys/dtypes.

            Returns:
                lab_dict (dict): Dict of lab frequencies/intensities (floats).
            """
            lab_dict = {}

            for _, datum in enumerate(data):
                # Assert valid key.
                key = datum.tag.split(namespace)[-1]
                if key not in lab_keys:
                    return None

                # Extract byte data to 32-bit floats.
                encoded_data = datum.text
                bytestring = base64.b64decode(encoded_data)
                lab_floats = array.array('f', bytestring)

                # Store in dictionary.
                lab_dict[key] = np.array(lab_floats)

            return lab_dict

        # Depending on what type of data, assign different keys/functions.
        if label == 'geometry':
            key_dict = geometry_keys
            func = parse_atom_to_dict

        elif label == 'transitions':
            key_dict = transition_keys
            func = parse_mode_to_dict

        elif label == 'laboratory':
            key_dict = laboratory_keys
            # No func here, since this XML tag is not repeated under
            # a single 'specie'. Handle it differently, as seen below.

        else:
            raise ValueError("No valid label, can't identify XML element.")

        # Isolate the XML tags for this data type.
        specie_tags = [x.tag for x in specie.getchildren()]
        try:
            relevant_indices = specie_tags.index(namespace + label)
        except ValueError:
            return None

        # Extract elements.
        data_block = specie.getchildren()[relevant_indices]
        data = data_block.getchildren()

        # Extract element subinformation and store in dictionary.
        if label == 'laboratory':
            lab_dictionary = parse_lab_to_dict(data, namespace, key_dict)
            return lab_dictionary

        else:
            container = []
            for _, datum in enumerate(data):
                data_dict = func(datum, namespace, key_dict)
                container.append(data_dict)

            return tuple(container)

    @classmethod
    def from_string(cls, string):
        """Generate XMLparser instance from a string."""
        # TODO: complete this method, if the general class structure is OK.
        return cls()
