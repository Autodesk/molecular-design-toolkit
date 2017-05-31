from __future__ import print_function, absolute_import, division
from future.builtins import *
from future import standard_library
standard_library.install_aliases()

# Copyright 2017 Autodesk Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import json
import zlib
import future.utils

from . import Alias

if future.utils.PY2:
    import dumbdbm
else:
    import dbm.dumb as dumbdbm


class CompressedJsonDbm(object):
    """ Quick-and-dirty interface to a DBM file
    """
    def __init__(self, filename, flag='r', dbm=dumbdbm):
        self.dbm = dbm
        if hasattr(dbm, 'open'):
            self.db = self.dbm.open(filename, flag)
        else:
            self.db = self.dbm(filename, flag)

    def __getattr__(self, item):
        return getattr(self.db, item)

    def __dir__(self):
        return list(self.__dict__.keys()) + dir(self.db)

    def __len__(self):
        return len(self.db)

    def __getitem__(self, key):
        gzvalue = self.db[key]
        return json.loads(zlib.decompress(gzvalue).decode())

    def __setitem__(self, key, value):
        gzvalue = zlib.compress(json.dumps(value))
        self.db[key] = gzvalue

    __contains__ = Alias('db.__contains__')


class ReadOnlyDumb(dumbdbm._Database):
    """ A read-only subclass of dumbdbm

    All possible operations that could result in a disk write have been turned into no-ops or raise
    exceptions
    """
    def _commit(self):
        # Does nothing!
        pass

    def __setitem__(self, key, value):
        raise NotImplementedError('This is a read-only database')

    def __delitem__(self, key):
        raise NotImplementedError('This is a read-only database')

    def _addkey(self, *args):
        assert False, 'Should never be here - this is a read-only database'


