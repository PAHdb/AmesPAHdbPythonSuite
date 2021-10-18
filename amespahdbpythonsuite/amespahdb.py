#!/usr/bin/env python3

import os
import copy
import hashlib
import tempfile
import time
from datetime import timedelta
import re

import pickle

from amespahdbpythonsuite.xmlparser import XMLparser


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

        if os.path.isfile("VERSION") and os.access("VERSION", os.R_OK):
            f = open("VERSION", 'r')
            dated = f.readline()
            f.close()
            self.message(f'SUITE DATED: {dated.upp}')

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

        if not os.path.isfile(filename) or not os.access(filename, os.R_OK):
            self.message(f'UNABLE TO REAS: {filename}')
            return

        self.__data = dict()

        self._joined = None

        self.parsefile(filename,
                       cache=keywords.get('cache', True),
                       check=keywords.get('check', True))

    def parsefile(self, filename, **keywords):
        """
        Method to parse or restore the database from cache.
        Called by :meth:`amespahdbpythonsuite.amespahdb.__init__`method.

        Parameters
        ----------
        filname : str
            String of the PAHdb filename and path.
        **keywords
            Arbitrary keyword arguments.

        """

        # Create MD5 hash function pickle.
        md5 = tempfile.gettempdir() + \
            '/' + hashlib.md5(open(filename, 'r').read().encode()).hexdigest()

        joined = md5 + '_joined.pkl'

        if not keywords.get('cache', True):
            if os.path.isfile(joined) and os.access(joined, os.R_OK):
                os.remove(joined)

        md5 += '.pkl'

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

    def pointer(self):
        return self.__data

    def __getkeybyuids(self, key, uids):
        """
        Get a dictionary of PAHdb properties
        retrieved by keyword for provided UIDs.

        Parameters
        ----------
        key : str
            Database keyword.
        uids : list of integers
            List of UIDs.

        Returns
        -------
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

        Parameters
        ----------
        uids : list of integers
            List of UIDs.

        Returns
        -------
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

        Parameters
        ----------
        uids : list of integers
            List of UIDs.

        Returns
        -------
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

        Parameters
        ----------
        uids : list of integers
            List of UIDs.

        Returns
        -------
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

        Parameters
        ----------
        uids : list of integers
            List of UIDs.

        Returns
        -------
        geometry object

        """
        from amespahdbpythonsuite.geometry import Geometry

        if type(uids) == int:
            uids = [uids]
        return Geometry(data=self.__getkeybyuids('geometry', uids), uids=uids)

    def search(self, query):
        """
        Search the database based on query input.

        Parameters
        ----------
        query : str
            String containing search query.

        Returns
        -------
        List of UIDs

        Example
        -------
        ``search('magnesium=0 oxygen=0 iron=0
        silicium=0 chx=0 ch2=0 c>20 h>0')``

        """

        return None

    def _tokenize(self, word):
        """
        A method called by :sec:`amespahdbpythonsuite.amespahdb.search`
        that creates a dictionary with keys based on the input word/category.

        Parameters
        ----------
        word : str

        """

        return None

    def _parsetokens(self, tokens):
        """
        Parse the dictionary of tokens created by
        :sec:`amespahdbpythonsuite.amespahdb.tokenize`
        and return string of expressions.

        Parameters
        ----------
        tokens : list
            List of dictionaries.

        Returns
        -------
        parsed : str
            String of expressions based on tokens.

        """

        return None

    def getversion(self):
        """
        Method to retrieve the PAHdb version.

        Returns
        -------
        String of PAHdb version.

        """
        return self.__data['version']

    def gettype(self):
        """
        Method to retrieve the PAHdb type.

        Returns
        -------
        String of PAHdb type.

        """
        return self.__data['type']

    def getdata(self):
        """
        Method to retrieve the database.

        Returns
        -------
        Dictionary containing the parsed database.

        """
        return self.__data

    @staticmethod
    def message(text, space=55):
        """
        A method to print terminal message.

        Parameters
        ----------
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
