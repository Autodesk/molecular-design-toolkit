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

import collections
from .descriptors import Alias

class Categorizer(dict):
    """
    Create a dict of lists from an iterable, with dict keys given by keyfn
    """

    def __init__(self, keyfn, iterable):
        super(Categorizer, self).__init__()
        self.keyfn = keyfn
        for item in iterable:
            self.add(item)

    def add(self, item):
        key = self.keyfn(item)
        if key not in self:
            self[key] = []
        self[key].append(item)


class ExclusiveList(object):
    """ Behaves like a list, but won't allow duplicate items with duplicate keys to be added.
    """
    def __init__(self, iterable=None, key=None):
        self._keys = collections.OrderedDict()
        if key is None:
            self._keyfn = self._identity
        else:
            self._keyfn = key

        if iterable is not None:
            self.extend(iterable)

    def append(self, obj):
        k = self._keyfn(obj)
        if k in self._keys:
            raise KeyError("'%s' can't be added because its key '%s' already exists" % (obj, k))
        else:
            self._keys[k] = obj

    def clear(self):
        self._keys = collections.OrderedDict()

    @staticmethod
    def _identity(obj):
        return obj

    def __iter__(self):
        return iter(self._keys.values())

    def __len__(self):
        return len(self._keys)

    def __getitem__(self, item):
        return list(self._keys.values())[item]

    def remove(self, obj):
        k = self._keyfn(obj)
        stored = self._keys[k]
        if obj is not stored:
            raise KeyError(obj)
        else:
            self._keys.pop(k)

    def extend(self, iterable):
        for item in iterable:
            self.append(item)

    def pop(self, index=None):
        if index is None:
            return self._keys.popitem()[1]
        else:
            k = list(self._keys.keys())[index]
            return self._keys.pop(k)

    def __repr__(self):
        return '%s(%s)' % (type(self).__name__, list(self._keys.values()))

    __str__ = __repr__


class DotDict(object):
    """ An attribute-accessible dictionary that preserved insertion order
    """

    def __delattr__(self, item):
        if not self.__dict__.get('_init', False):
            super().__delattr__(item)
        else:
            try:
                del self._od[item]
            except KeyError:
                raise AttributeError()

    def __delitem__(self, key):
        if not self.__dict__.get('_init', False):
            super().__delitem__(key)
        else:
            del self._od[key]

    def __reduce__(self):
        return DotDict, self._od

    def copy(self):
        return self.__class__(self._od.copy())

    def __init__(self, *args, **kwargs):
        self._od = collections.OrderedDict(*args, **kwargs)
        self._init = True

    def __eq__(self, other):
        try:
            return self._od == other._od
        except AttributeError:
            return False

    def __repr__(self):
        return str(self._od).replace('OrderedDict', self.__class__.__name__)

    def __getattr__(self, key):
        if not self.__dict__.get('_init', False):
            return self.__getattribute__(key)
        if key in self._od:
            return self._od[key]
        else:
            raise AttributeError(key)

    def __setattr__(self, key, val):
        if not self.__dict__.get('_init', False):
            super().__setattr__(key, val)
        else:
            self._od[key] = val

for _v in ('keys values items __iter__ __getitem__  __len__ __contains__ '
           '__setitem__ pop setdefault get update').split():
    setattr(DotDict, _v, Alias('_od.%s' % _v))
