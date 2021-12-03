#!/usr/bin/env python3

from amespahdbpythonsuite.amespahdb import AmesPAHdb
message = AmesPAHdb.message


class Data(object):
    """
    AmesPAHdbPythonSuite data class

    """

    def __init__(self, d=None, **keywords):
        """
        Initialize Data class.

        """

        self.set(d, **keywords)

    def set(self, d=None, **keywords):
        """
        Populate data dictionary.

        """
        if d:
            # Check if expected keywords are present in provided dictionary,
            # otherwise assign them to instance variables.
            if d.get('type', '') == self.__class__.__name__:
                if not keywords.get('type'):
                    self.type = d['database']
                if not keywords.get('version'):
                    self.version = d['version']
                if not keywords.get('data'):
                    self.data = d['data']
                if not keywords.get('uids'):
                    self.uids = d['uids']
                    print('no uids keyword')
                if not keywords.get('model'):
                    self.model = d['model']
                if not keywords.get('units'):
                    self.units = d['units']

        # Match keywords of provided dictionary to corresponding instance variables.
        if keywords.get('type'):
            self.type = keywords.get('type')
        if keywords.get('version'):
            self.version = keywords.get('version')
        if keywords.get('data'):
            self.data = keywords.get('data')
        if keywords.get('pahdb'):
            self.pahdb = keywords.get('pahdb')
        if keywords.get('uids'):
            self.uids = keywords.get('uids')
        if keywords.get('model'):
            self.model = keywords.get('model')
        if keywords.get('units'):
            self.units = keywords.get('units')

        # Check for database and versioning mismatch between provided dictionary and parsed database.
        if self.pahdb:
            if self.pahdb['database'] != self.type:

                message(f'DATABASE MISMATCH: {self.pahdb["database"]} != {self.type}')
                return

            if self.pahdb['version'] != self.version:

                message(f'VERSION MISMATCH: {self.pahdb["version"]} != {self.version}')
                return

    def get(self):
        """
        Return data dictionary with expected keywords.

        """
        return {'type': self.__class__.__name__,
                'database': self.type,
                'version': self.version,
                'data': self.data,
                'uids': self.uids,
                'model': self.model,
                'units': self.units}

    def intersect(self, uids):
        """
        Updates data to the intersection with provided UIDs.

        Parameters
        ----------
        uids : list of integers

        """
        self.uids = set(self.uids) & set(uids)

        count = len(self.uids)

        if count == 0:

            message('NO INTERSECTION FOUND')

            return

        message(f'INTERSECTION FOUND: {count}')

        self.data = {key: self.data[key] for key in self.uids}
