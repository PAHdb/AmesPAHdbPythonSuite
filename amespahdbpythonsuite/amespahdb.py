#!/usr/bin/env python3

import os
import sys
import copy
import hashlib
import tempfile
import time
import re
from datetime import timedelta

import pickle

from amespahdbpythonsuite.xmlparser import XMLparser
import amespahdbpythonsuite as suite


class AmesPAHdb:
    """
    AmesPAHdbPythonSuite main class.
    Contains methods to parse the database, perform search based on query,
    and retrieve UIDs.
    Calls methods:
    :meth:`amespahdbpythonsuite.Transitions`,
    :meth:`amespahdbpythonsuite.Laboratory`,
    :meth:`amespahdbpythonsuite.Species`,
    :meth:`amespahdbpythonsuite.Geometry`,
    to retrieve the respective objects.

    """

    def __init__(self, **keywords):
        """
        Initialize amespahdbpythonsuite class. Prints basic PAHdb info
        and calls the :meth:`amespahdbpythonsuite.amespahdb.parsefile`
        method to parse or restore the database.

        """

        intro = ['AmesPAHdbPythonSuite', 'by', 'Dr. Christiaan Boersma', 'and',
                 'Dr. Alexandros Maragkoudakis']
        self.message(intro)

        self.message(f'SUITE VERSION: {suite.__version__}')

        self.message('WEBSITE: HTTP://WWW.ASTROCHEM.ORG/PAHDB/')
        self.message('CONTACT: CHRISTIAAN.BOERSMA@NASA.GOV')

        filename = keywords.get('filename')
        if not filename:
            filename = os.environ.get("AMESPAHDEFAULTDB")
            if not filename:
                msg = ['DATABASE NOT FOUND:',
                       'SET SYSTEM AMESPAHDEFAULTDB',
                       'ENVIRONMENT VARIABLE']
                self.message(' '.join(msg))
                # TODO: Turn the sys.exit into exceptions.
                sys.exit(1)

        if not filename or not os.path.isfile(filename) or not os.access(filename, os.R_OK):
            self.message(f'UNABLE TO READ: {filename}')
            # TODO: Turn the sys.exit into exceptions.
            sys.exit(2)

        self.__data = {}

        self._joined = None

        self.parsefile(filename,
                       cache=keywords.get('cache', True),
                       check=keywords.get('check', True))

    def parsefile(self, filename, **keywords):
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
        md5 = tempfile.gettempdir() + \
            '/' + hashlib.md5(open(filename, 'r').read().encode()
                              ).hexdigest() + '.pkl'

        # Check if database is dumped in cache and restore it.
        if (keywords.get('cache', True) and os.path.isfile(md5) and os.access(md5, os.R_OK)):

            self.message('RESTORING DATABASE FROM CACHE')

            # Start timer
            tstart = time.perf_counter()

            # Open and read the dumped database.
            with open(md5, 'rb') as f:
                self.__data = pickle.load(f)

            # Store the dumped database filename.
            self.__data['filename'] = md5

            # stop timer and calculate elapsed time
            elapsed = timedelta(seconds=(time.perf_counter() - tstart))

            info = [f'FILENAME                    : {md5}',
                    f'ORIGNINAL FILENAME          : {filename}',
                    f'PARSE TIME                  : {elapsed}',
                    f'VERSION (DATE)              : {self.__data["version"]} ({self.__data["date"]})',
                    f'COMMENT                     : {self.__data["comment"]}']
            self.message(info, space=0)

        else:

            self.message('PARSING DATABASE: THIS MAY TAKE A FEW MINUTES')

            # Start timer
            tstart = time.perf_counter()

            # Call XMLparser module to parse the database.
            parser = XMLparser(filename=filename,
                               validate=keywords.get('check', True))

            # Store the database into self.__self.data.
            self.__data = parser.to_pahdb_dict()

            # Dump the database into pickle in the cache directory.
            with open(md5, 'wb') as f:
                pickle.dump(self.__data, f, pickle.HIGHEST_PROTOCOL)

            # Store the dumped database filename.
            self.__data['filename'] = md5

            # stop timer and calculate elapsed time
            elapsed = timedelta(seconds=(time.perf_counter() - tstart))
            print(f'Elapsed time: {elapsed}\n')

    def __getkeybyuids(self, key, uids):
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
        # TODO Handle empty lists or invalid UIDs.
        if key == 'species':
            return copy.deepcopy(dict((uid, self.__data['species'][uid])
                                      for uid in uids if uid in self.__data['species'].keys()))
        else:
            return copy.deepcopy(dict((uid, self.__data['species'][uid][key])
                                      for uid in uids if uid in self.__data['species'].keys()))

    def gettransitionsbyuid(self, uids):
        """
        Retrieve and return transitions object based on UIDs input.
        UIDs should be a list, e.g. the output of search method.
        Calls the :meth:`amespahdbpythonsuite.transitions.Transitions` class.

        Parameters:
            uids : list of integers
                List of UIDs.

        Returns:
            transitions object

        """
        from amespahdbpythonsuite.transitions import Transitions
        if type(uids) == int:
            uids = [uids]

        return Transitions(type=self.__data['database'],
                           version=self.__data['version'],
                           data=self.__getkeybyuids('transitions', uids),
                           pahdb=self.__data,
                           uids=uids,
                           model={'type': 'zerokelvin_m',
                                  'temperature': 0.0,
                                  'description': ''},
                           units={'abscissa': {'unit': 1,
                                               'str': 'frequency [wavenumber]'},
                                  'ordinate': {'unit': 2,
                                               'str': 'integrated cross-section' + '[km/mol]'}}
                           )

    def getlaboratorybyuid(self, uids):
        """
        Retrieve and return laboratory database object based on UIDs input.
        UIDs should be a list, e.g. the output of search method.
        Calls the :meth:`amespahdbpythonsuite.laboratory.Laboratory` class.

        Parameters:
            uids : list of integers
                List of UIDs.

        Returns:
            laboratory database object

        """
        # Check if the experimental database is loaded.
        if self.__data['database'] != 'experimental':

            self.message('EXPERIMENTAL DATABASE REQUIRED')

            return None

        from amespahdbpythonsuite.laboratory import Laboratory

        if type(uids) == int:
            uids = [uids]

        return Laboratory(type=self.__data['database'],
                          version=self.__data['version'],
                          data=self.__getkeybyuids('laboratory', uids),
                          pahdb=self.__data,
                          uids=uids,
                          model={'type': 'laboratory_m',
                                 'temperature': 0.0,
                                 'description': ''},
                          units={'abscissa': {'unit': 1,
                                              'str': 'frequency [wavenumber]'},
                                 'ordinate': {'unit': 2,
                                              'str': 'absorbance' + '[-log(I/I$_{0})$'}}
                          )

    def getspeciesbyuid(self, uids):
        """
        Retrieve and return species object based on UIDs input.
        UIDs should be a list, e.g. the output of search method.
        Calls the :meth:`amespahdbpythonsuite.species.Species` class.

        Parameters:
            uids : list of integers
                List of UIDs.

        Returns:
            species object

        """
        # TODO: handle UIDs without 'references'.
        from amespahdbpythonsuite.species import Species

        if type(uids) == int:
            uids = [uids]

        return Species(type=self.__data['database'],
                       version=self.__data['version'],
                       data=self.__getkeybyuids('species', uids),
                       pahdb=self.__data,
                       uids=uids,
                       references=self.__getkeybyuids('references', uids),
                       comments=self.__getkeybyuids('comments', uids)
                       )

    def getgeometrybyuid(self, uids):
        """
        Retrieve and return geometry object based on UIDs input.
        UIDs should be a list, e.g. the output of search method.
        Calls the :meth:`amespahdbpythonsuite.geometry.Geometry` class
        and :meth:`amespahdbpythonsuite.amespahdb.__getkeybyuids` method.

        Parameters:
            uids : list of integers
                List of UIDs.

        Returns:
            geometry object

        """
        from amespahdbpythonsuite.geometry import Geometry

        if type(uids) == int:
            uids = [uids]
        return Geometry(data=self.__getkeybyuids('geometry', uids), uids=uids)

    def search(self, query):
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
            return

        words = []

        n = len(query)

        i = 0
        while True:

            if i == n:
                break

            while query[i] == ' ':
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
                while i < n and query[i] not in [" ", "=", "<", ">", "&", "|", "(", ")", "!"]:
                    token += query[i]
                    i += 1

            words.append(token)

        tokens = []

        for word in words:
            tokens.append(self._tokenize(word))

        code = self._parsetokens(tokens)

        if not code:
            return None

        found = eval(
            f"[item[0] for item in _AmesPAHdb__data['species'].items() if ({code})]", self.__dict__)

        return found

    def _tokenize(self, word):
        """
        A method called by :sec:`amespahdbpythonsuite.amespahdb.search`
        that creates a dictionary with keys based on the input word/category.

        Parameters:
            word : str

        """

        word = word.lower()

        token = {}

        charge = {'anion': 'item[1]["charge"] < 0',
                  'cation': 'item[1]["charge"] > 0',
                  'neutral': 'item[1]["charge"] == 0',
                  'positive': 'item[1]["charge"] > 0',
                  'negative': 'item[1]["charge"] < 0',
                  '-': 'item[1]["charge"] == -1',
                  '+': 'item[1]["charge"] == 1',
                  '++': 'item[1]["charge"] == 2',
                  '+++': 'item[1]["charge"] == 3',
                  '---': 'item[1]["charge"] == -3'}

        identities = {'uid': 'item[0]',
                      'identifier': 'item[0]',
                      'atoms': 'len(item[1]["geometry"])',
                      'hydrogen': 'len([c for c in item[1]["geometry"] if c["type"] == 1])',
                      'carbon': 'len([c for c in item[1]["geometry"] if c["type"] == 6])',
                      'nitrogen': 'len([c for c in item[1]["geometry"] if c["type"] == 7])',
                      'oxygen': 'len([c for c in item[1]["geometry"] if c["type"] == 8])',
                      'magnesium': 'len([c for c in item[1]["geometry"] if c["type"] == 12])',
                      'silicium': 'len([c for c in item[1]["geometry"] if c["type"] == 14])',
                      'iron': 'len([c for c in item[1]["geometry"] if c["type"] == 26])',
                      'h': 'len([c for c in item[1]["geometry"] if c["type"] == 1])',
                      'c': 'len([c for c in item[1]["geometry"] if c["type"] == 6])',
                      'n': 'len([c for c in item[1]["geometry"] if c["type"] == 7])',
                      'o': 'len([c for c in item[1]["geometry"] if c["type"] == 8])',
                      'mg': 'len([c for c in item[1]["geometry"] if c["type"] == 12])',
                      'si': 'len([c for c in item[1]["geometry"] if c["type"] == 14])',
                      'fe': 'len([c for c in item[1]["geometry"] if c["type"] == 26])',
                      # TODO make transition search work
                      # 'wavenumber': 'item[1]["transitions"]["frequency"]',
                      # 'absorbance': 'item[1]["transitions"]["intensity"]',
                      # 'frequency': 'item[1]["transitions"]["frequency"]',
                      # 'intensity': 'item[1]["transitions"]["intensity"]',
                      'ch2': 'item[1]["n_ch2"]',
                      'chx': 'item[1]["n_chx"]',
                      'solo': 'item[1]["n_solo"]',
                      'duo': 'item[1]["n_duo"]',
                      'trio': 'item[1]["n_trio"]',
                      'quartet': 'item[1]["n_quartet"]',
                      'quintet': 'item[1]["n_quintet"]',
                      'charge': 'item[1]["charge"]',
                      'symmetry': 'item[1]["symmetry"]',
                      'weight': 'item[1]["weight"]',
                      'scale': 'item[1]["scale"]',
                      'energy': 'item[1]["total_e"]',
                      'zeropoint': 'item[1]m["vib_e"]',
                      'experiment': 'item[1]["exp"]'}

        logical = {'and': 'and', 'or': 'or', '|': 'or', '&': 'and'}

        comparison = {'<': '<', 'lt': '<', '>': '>', 'gt': '>', '=': '==', 'eq': '==',
                      '<=': '<=', 'le': '<=', '>=': '>=', 'ge': '>=', 'with': 'and', 'ne': '!=', '!=': '!='}

        transfer = {'(': '(', ')': ')'}

        if word.isnumeric():
            token['type'] = "NUMERIC"
            token['translation'] = word
        elif word in charge:
            token['type'] = "CHARGE"
            token['translation'] = charge[word]
        elif word in identities:
            token['type'] = "IDENTITY"
            token['translation'] = identities[word]
        elif word in logical:
            token['type'] = "LOGICAL"
            token['translation'] = logical[word]
        elif word in comparison:
            token['type'] = "COMPARISON"
            token['translation'] = comparison[word]
        elif word in transfer:
            token['type'] = "TRANSFER"
            token['translation'] = transfer[word]
        elif re.search(
            '(mg+|si+|fe+|[chno]+)([0-9]*)(mg+|si+|fe+|[chno]+)([0-9]*)(mg+|si+|fe+|[chno]*)([0-9]*)',
                word):
            token['type'] = "FORMULA"
            token['translation'] = f"item[1]['formula'] == '{word.upper()}'"
        else:
            # TODO: add search by compound name
            token['type'] = "IGNORE"
            token['translation'] = word

        token['valid'] = True

        return token

    def _parsetokens(self, tokens):
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

        prev = None

        current = 0

        if ntokens > 1:
            next = 1
        else:
            next = None

        parsed = ''

        while current is not None:

            if tokens[current]['type'] == 'FORMULA':
                if prev is not None:
                    if not (tokens[prev]['type'] != 'LOGICAL' and tokens[prev]['valid']):
                        parsed += ' or '
                parsed += tokens[current]['translation']
            elif tokens[current]['type'] == 'IDENTITY':
                if prev is not None:
                    if not (tokens[prev]['type'] == 'LOGICAL' or
                       tokens[prev]['type'] == 'TRANSFER' or
                       tokens[prev]['type'] == 'TRANSFER' and tokens[prev]['valid']):
                        parsed += ' and '
                if next is not None:
                    if tokens[next]['type'] == 'COMPARISON':
                        parsed += ' ' + tokens[current]['translation']
                    else:
                        parsed += ' ' + tokens[current]['translation'] + ' > 0'
            elif tokens[current]['type'] == 'NUMERIC':
                if prev is not None:
                    if tokens[prev]['type'] == 'COMPARISON' and tokens[prev]['valid']:
                        parsed += ' ' + tokens[current]['translation']
                    else:
                        tokens[current].valid = False
            elif tokens[current]['type'] == 'LOGICAL':
                if prev is not None:
                    if (tokens[prev]['type'] == 'IDENTITY' or tokens[prev]['type'] == 'NUMERIC' or
                        tokens[prev]['type'] == 'FORMULA' or tokens[prev]['type'] == 'CHARGE' and
                            tokens[prev]['valid']):
                        if next is not None:
                            if (tokens[next]['type'] == 'TRANSFER'):
                                parsed += tokens[current]['translation']
                            elif (tokens[next]['type'] == 'IDENTITY' or
                                  tokens[next]['type'] == 'NUMERIC' or
                                  tokens[next]['type'] == 'FORMULA' or
                                  tokens[next]['type'] == 'CHARGE'):
                                parsed += ' ' + tokens[current]['translation']
                            else:
                                tokens[current].valid = False
            elif tokens[current]['type'] == 'COMPARISON':
                if prev is not None:
                    if tokens[prev]['type'] == 'IDENTITY' and tokens[prev]['valid']:
                        if next is not None:
                            if tokens[next]['type'] == 'NUMERIC':
                                parsed += ' ' + tokens[current]['translation']
                            else:
                                tokens[current]['valid'] = False
            elif tokens[current]['type'] == 'CHARGE':
                if prev is not None:
                    if not (tokens[prev]['type'] == 'LOGICAL' and tokens[prev]['valid']):
                        parsed += ' and '
                parsed += ' ' + tokens[current]['translation']
            elif tokens[current]['type'] == 'TRANSFER':
                parsed += tokens[current]['translation']
            elif tokens[current]['type'] == 'NAME':
                if prev is not None:
                    if not (tokens[prev]['type'] == 'LOGICAL' and tokens[prev]['valid']):
                        parsed += ' and '
                # TODO implement name
                # parsed += f"item[1]['comments'] == {tokens[current]['translation']}"
            elif tokens[current]['type'] == 'IGNORE':
                print(f"'{tokens[current]['translation']}' NOT UNDERSTOOD")

                return ''

            prev = current

            current = next

            if next:
                if next == ntokens - 1:
                    next = None
                else:
                    next += 1

        return parsed

    def getversion(self):
        """
        Method to retrieve the PAHdb version.

        Returns:
            String of PAHdb version.

        """
        return self.__data['version']

    def gettype(self):
        """
        Method to retrieve the PAHdb type.

        Returns:
            String of PAHdb type.

        """
        return self.__data['type']

    def getdata(self):
        """
        Method to retrieve the database.

        Returns:
            Dictionary containing the parsed database.

        """
        return self.__data

    @staticmethod
    def message(text, space=55):
        """
        A method to print terminal message.

        Parameters:
            text : string or list of strings.
                Text to be displayed.
            space : integer
                Number to indent the text.

        """
        line = 57 * '='
        print(line)
        if type(text) is list:
            for t in text:
                print(t.center(space))
        else:
            print(text.center(space))
        print(line)
