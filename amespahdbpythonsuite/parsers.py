#!/usr/bin/env python3
"""
parsers.py

Parse a PAHdb XML file, with or without schema checking.
"""

import array
import base64
import os
import urllib.request

from lxml import etree as ElementTree
import numpy as np


class XMLparser:
    """Parse a PAHdb XML file.

    Optional behaviour includes validating the XML schema.

    Attributes:
        check_schema: Whether to validate the XML schema or not.
        filename: XML file location. If not given, it is assigned to the
            "AMESPAHDEFAULTDB" environment variable.

    Note:
        `validated` (bool) indicates if the schema is valid (True) or
            not (False). Defaults to False until validation is performed.
    """

    def __init__(self, filename=None, check_schema=True):
        """Inits XMLparser with schema checking on, no given filename."""
        self.filename = filename
        self.check_schema = check_schema
        self.validated = False  # Assume schema has not yet been checked.

        # Identify location of PAHdb XML file.
        if self.filename is None:
            try:
                self.filename = os.environ.get("AMESPAHDEFAULTDB")
            except OSError as e:
                print("Must set AMESPAHDEFAULTDB environment variable.")
                raise e

        # Parse XML file to dictionary.
        self._database, self._info = self._parse()

    def __repr__(self):
        """Class representation."""
        return (f'{self.__class__.__name__}('
                f'check_schema={self.check_schema!r}, '
                f'filename={self.filename!r}, '
                f'validated={self.validated!r})')

    def __str__(self):
        """A description of the object."""
        return (f'An XML parser from the Ames PAHdb Python Suite. '
                f'XML file: \n{self.filename!r}.')

    def _parse(self):
        """Perform the actual parsing, with or without validation.

        Returns:
            database (dict): Dictionary, with the UIDs as keys,
                containing the transitions, geometry data, as well as
                UID metadata, references, comments.
            info (dict): Dictionary containing general information about
                the XML file, including filename, type, date, full
                (bool), version, comment, and number of species.
        """

        if self.check_schema:
            tree = ElementTree.parse(self.filename)
            root = tree.getroot()
            self._validate_schema(tree, root)
            formatted_tree = \
                ElementTree.iterwalk(tree, events=("start", "end"))

        else:
            formatted_tree = \
                ElementTree.iterparse(self.filename, events=("start", "end"))

        database, info = self._tree_to_dict(formatted_tree)

        return database, info

    def _validate_schema(self, tree, root):
        """Validate XML schema.

        Args:
            tree (ElementTree.ElementTree): parsed XML file.
            root: root element of ElementTree tree.

        Returns:
            N/A. If the XML schema is successfully validated, the instance
            variable self.validated is set to True (bool).
        """

        schema = root.get('{http://www.w3.org/2001/XMLSchema-instance}' +
                          'schemaLocation')
        if schema:
            _, uri = schema.split(' ', 1)
            doc = ElementTree.parse(urllib.request.urlopen(uri))
            xmlschema = ElementTree.XMLSchema(doc)
            try:
                xmlschema.assertValid(tree)
            except Exception as e:
                raise e
            else:
                self.validated = True

    def _tree_to_dict(self, tree):
        """Convert the tree to a dict of dicts.

        Args:
            tree: lxml tree of PAHdb file.

        Returns:
            database_dict: Dictionary, with the UIDs as keys,
                containing the transitions, geometry data, as well as
                UID metadata, references, comments.
            info_dict: Dictionary containing general information about
                the XML file, including filename, type, date, full
                (bool), version, comment, and number of species.
        """
        # Advance to root and identify namespace.
        event, root = next(tree)
        namespace = self.determine_xml_namespace(root)

        # Root dict attributes.
        root_dict = root.attrib

        # Form important dicts, to be returned.
        database_dict = {}
        info_dict = {}
        first_run = True

        # Iterate through the tree, clearing elements as we go.
        for _, (event, elem) in enumerate(tree):

            # Store the top XML comment.
            if first_run:
                if event == "end" and elem.tag == namespace + 'comment':
                    pahdb_comment = elem.text
                    first_run = False

            # Extract UID data on a per-tag ('specie') basis.
            if event == "end" and elem.tag == namespace + 'specie':
                # Extract soft data.
                uid, metadata, comments, reference = \
                    self._extract_info(elem, namespace)

                # Extract hard data.
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

    def to_dict(self):
        """Return the database and info dictionaries."""
        return self._database, self._info

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

        def get_metadata_keys():
            """XML tags for the children of a <specie> element,
            metadata only."""
            metadata_keys = {
                'references': str,
                'comments': str,
                'formula': str,
                'charge': int,
                'symmetry': str,
                'weight': float,
                'total_e': float,
                'vib_e': float,
                'method': str,
                'n_solo': int,
                'n_duo': int,
                'n_trio': int,
                'n_quartet': int,
                'n_quintet': int,
                'n_ch2': int,
                'n_chx': int,
            }
            return metadata_keys

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

            # Retrieve metadata from each child.
            metadata_keys = get_metadata_keys()

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

        def geometry_keys():
            """Returns the XML tags for the children of <geometry> as a
            dictionary."""
            geo_keys = {
                'position': int,
                # 'uid': int,
                # 'type': int,
                'x': float,
                'y': float,
                'z': float,
            }
            return geo_keys

        def transition_keys():
            """Returns the XML tags for the children of <transitions> as a
            dictionary."""
            tran_keys = {
                'frequency': float,
                # 'uid': int,
                'scale': float,
                'intensity': float,
                'symmetry': str,
            }
            return tran_keys

        def laboratory_keys():
            """Returns the XML tags for the children of <laboratory> as a
            dictionary."""
            lab_keys = {
                'frequency': float,
                'intensity': float,
            }
            return lab_keys

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
            key_dict = geometry_keys()
            func = parse_atom_to_dict

        elif label == 'transitions':
            key_dict = transition_keys()
            func = parse_mode_to_dict

        elif label == 'laboratory':
            key_dict = laboratory_keys()
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
