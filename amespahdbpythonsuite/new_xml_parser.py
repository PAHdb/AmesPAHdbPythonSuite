#!/usr/bin/env python3
"""
new_xml_parser.py

Try a new parser.
"""

from pkg_resources import resource_filename

from ipdb import set_trace as st
from lxml import etree


def species_handler(action, elem, context):
    """Parse a PAHdb XML <species> tag."""
    database = {}

    while True:
        action, elem = next(context)
        tag = etree.QName(elem).localname

        if action == 'start' and tag == 'specie':
            uid, specie_dict = specie_handler(action, elem, context)
            database[uid] = specie_dict

        elif action == 'end' and tag == 'species':
            break

        del action, elem

    return database


def specie_handler(action, elem, context):
    """Parse a PAHdb XML <specie> tag."""

    def specie_geometry_handler(elem, context):
        """<specie> tag: Parse its child <geometry> tag."""
        geometry = []

        while True:
            action, elem = next(context)
            tag = etree.QName(elem).localname

            if action == 'end' and tag == 'geometry':
                break

            if action != 'start':
                pass

            if tag == 'atom':
                atom_dict = {}
                while True:
                    action, elem = next(context)
                    tag = etree.QName(elem).localname
                    if action == 'end' and tag == 'atom':
                        print("broke out")
                        break
                    print(elem.text, elem.tag, elem.attrib)
                    atom_dict[tag] = float(elem.text)
                    del action, elem

                geometry.append(atom_dict)

            del action, elem

        return geometry

    def specie_transitions_handler(elem, context):
        """<specie> tag: Parse its child <transitions> tag."""
        transitions = []

        while True:
            action, elem = next(context)
            tag = etree.QName(elem).localname

            if action == 'end' and tag == 'transitions':
                break

            if action != 'start':
                pass

            if tag == 'mode':
                mode_dict = {}
                while True:
                    action, elem = next(context)
                    tag = etree.QName(elem).localname
                    if action == 'end' and tag == 'mode':
                        break

                    if elem.attrib:
                        mode_dict.update(elem.attrib)

                    try:
                        # TODO: fix this thang
                        # print(elem.text, elem.tag)
                        # if not elem.text:
                        #     st()
                        value = float(elem.text)
                    except ValueError:
                        value = elem.text

                    mode_dict[tag] = value

                    del action, elem

                transitions.append(mode_dict)

            del action, elem

        return transitions

    uid = int(elem.attrib['uid'])
    specie_dict = {}
    metadata_dict = {}

    while True:

        action, elem = next(context)
        tag = etree.QName(elem).localname

        if action == 'end' and tag == 'specie':
            break

        if action != 'start':
            pass

        # Parse comments
        if tag == 'comments':
            comments = []
            while True:
                action, elem = next(context)
                tag = etree.QName(elem).localname
                if action == 'start' and tag == 'comment':
                    comments.append(elem.text)
                elif action == 'end' and tag == 'comments':
                    break
                del action, elem

            specie_dict['comments'] = tuple(comments)

        # Parse references
        elif tag == 'references':
            references = []
            while True:
                action, elem = next(context)
                tag = etree.QName(elem).localname
                if action == 'start' and tag == 'reference':
                    comments.append(elem.text)
                elif action == 'end' and tag == 'references':
                    break
                del action, elem

            specie_dict['references'] = tuple(references)

        # Parse geometry
        elif tag == 'geometry':
            specie_dict['geometry'] = specie_geometry_handler(elem, context)

        # Parse transitions
        elif tag == 'transitions':
            specie_dict['transitions'] = \
                specie_transitions_handler(elem, context)

        else:
            try:
                value = float(elem.text)
                print(value, elem.tag, elem.attrib, elem.text)
                if value % 1 == 0:
                    value = int(value)
            except ValueError:
                value = elem.text
            metadata_dict[tag] = value

        del action, elem

    specie_dict['metadata'] = metadata_dict

    return uid, specie_dict





sample_xml_file = resource_filename(
    # 'amespahdbpythonsuite', 'resources/pahdb-theoretical_cutdown.xml'
    'amespahdbpythonsuite', 'resources/pahdb-complete-theoretical-v3.10.xml'
)

context = etree.iterparse(sample_xml_file,
                          events=("start", "end"))


info = {}

while True:
    action, elem = next(context)
    print(action, elem)

    if action == 'end' and tag == 'pahdatabase':
        break

    tag = etree.QName(elem).localname
    print("%s: %s; %s, %s" % (action, elem.tag, elem.attrib, elem.text))

    if action == 'start' and tag == 'pahdatabase':
        info.update(elem.attrib)

    elif action == 'start' and tag == 'comment':
        info['comment'] = elem.text

    elif action == 'start' and tag == 'species':
        database = species_handler(action, elem, context)

    del action, elem

info['nspecies'] = len(database.keys())



# # desired output:
# info_dict = {
#     'filename': full_path,
#     'type': root_dict['database'],
#     'date': root_dict['date'],
#     'full': bool(root_dict['full']),
#     'version': root_dict['version'],
#     'comment': pahdb_comment,
#     'nspecies': len(db_keys),
# }

# database_dict[uid] = {
#     'metadata': metadata,
#     'geometry': geometry_dict,
#     'transitions': transitions_dict,
#     'comments': comments,
#     'reference': reference,
#     'laboratory': laboratory_dict,
# }


