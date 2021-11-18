#!/usr/bin/env python3

from amespahdbpythonsuite.amespahdb import AmesPAHdb
message = AmesPAHdb.message


class Data():
    """
    AmesPAHdbPythonSuite data class

    """

    def __init__(self, d=None, **keywords):
        self.set(d, **keywords)

    def set(self, d=None, **keywords):
        if d:
            if d.get('type', '') == self.__class__.__name__:
                if not keywords.get('type'):
                    self.type = d['database']
                if not keywords.get('version'):
                    self.version = d['version']
                if not keywords.get('data'):
                    self.data = d['data']
                if not keywords.get('uids'):
                    self.uids = d['uids']
                if not keywords.get('model'):
                    self.model = d['model']
                if not keywords.get('units'):
                    self.units = d['units']

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

        if self.pahdb:

            if self.pahdb['database'] != self.type:

                message(f'DATABASE MISMATCH: {self.pahdb["type"]} != {self.type}')
                return

            if self.pahdb['version'] != self.version:

                message(f'VERSION MISMATCH: {self.pahdb["version"]} != {self.version}')
                return

    def get(self):
        return {'type': self.__class__.__name__,
                'database': self.type,
                'version': self.version,
                'data': self.data,
                'uids': self.uids,
                'model': self.model,
                'units': self.units}
