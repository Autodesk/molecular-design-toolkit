# Copyright 2016 Autodesk Inc.
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
from bisect import bisect_left, bisect_right
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


class AutoIndexList(list):
    """ Every object in the list has an automatically assigned `index` attribute.
    """
    def __init__(self, *args, **kwargs):
        super(AutoIndexList, self).__init__(*args, **kwargs)
        self._reindex()

    def _reindex(self):
        for idx, item in enumerate(self):
            item.index = idx

    def append(self, o):
        o.index = len(self)
        super(AutoIndexList, self).append(o)

    def extend(self, l):
        for i, o in enumerate(l):
            o.index = len(self) + i
        super(AutoIndexList, self).extend(l)

    def pop(self, *args):
        result = super(AutoIndexList, self).pop(*args)
        if args:
            self._reindex()
        return result


# Always reindex after using any of these methods:
def make_wrapper(name):
    def w(self, *args, **kwargs):
        result = getattr(list, name)(self, *args, **kwargs)
        self._reindex()
        return result
    w.__doc__ = getattr(list, name).__doc__
    w.__name__ = name
    return w

for _wrapped in ['__setitem__', '__setslice__', 'remove', 'sort', 'insert']:
    setattr(AutoIndexList, _wrapped, make_wrapper(_wrapped))


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
        return self._keys.itervalues()

    def __len__(self):
        return len(self._keys)

    def __getitem__(self, item):
        result = self._keys.values()[item]
        if isinstance(item, slice):
            return self.__class__(result, self._keyfn)
        else:
            return result

    def __setitem__(self, item, value):
        raise NotImplementedError()

    def insert(self, index, item):
        raise NotImplementedError()

    def remove(self, obj):
        k = self._keyfn(obj)
        stored = self._keys[k]
        if obj is not stored:
            raise ValueError(obj)
        else:
            self._keys.pop(k)

    def extend(self, iterable):
        for item in iterable:
            self.append(item)

    def pop(self, index=None):
        if index is None:
            return self._keys.popitem()[1]
        else:
            k = self._keys.keys()[index]
            return self._keys.pop(k)

    def __repr__(self):
        return '%s(%s)' % (type(self).__name__, self._keys.values())

    __str__ = __repr__


class DotDict(object):
    """Class for use as a dictionary with attribute-style access.

    This class is our explicit container for storing data in an object's __dict__.
    For convenience, it also supports item-style access.

    These classes are provided as an IDE feature - real-time autocomplete in ipython is much
    easier to use with attribute access. They little-to-no value in other environments, but
    are not actively harmful.

    All the standard dictionary functions can be accessed at this object's __dict__ attribute.
    """
    def __init__(self, *args, **kwargs):
        self.update(dict(*args, **kwargs))

    def __getattr__(self, item):
        return getattr(self.__dict__, item)

    def __repr__(self):
        return '%s(%s)' % (self.__class__.__name__, repr(self.__dict__))

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def copy(self):
        return self.__class__(self.iteritems())

# add dictionary magic method aliases
for _delegate in ['__getitem__', '__setitem__', '__contains__', '__iter__',
                  '__delitem__', '__len__']:
    setattr(DotDict, _delegate, Alias('__dict__.%s' % _delegate))


# TODO: just subclass OrderedDict, come on now
class OrderedDotDict(DotDict):
    """Dict with items accessible as attributes"""
    def __init__(self, *args, **kwargs):
        self.__keyorder__ = collections.OrderedDict()
        self.update(*args, **kwargs)

    __iter__ = Alias('__keyorder__.__iter__')
    keys = Alias('__keyorder__.keys')
    __len__ = Alias('__keyorder__.__len__')

    iterkeys = __iter__

    def itervalues(self):
        for k in self:
            yield self[k]

    def iteritems(self):
        for k in self:
            yield k, self[k]

    def __repr__(self):
        return '%s({%s})'%(self.__class__.__name__, self.items())

    def values(self):
        return list(self.itervalues())

    def items(self):
        return list(self.iteritems())

    def __setitem__(self, key, value):
        super(OrderedDotDict, self).__setitem__(key, value)
        self.__keyorder__[key] = None

    def __setattr__(self, key, value):
        super(OrderedDotDict, self).__setattr__(key, value)
        if key[:2] != '__':
            self.__keyorder__[key] = None

    def __delitem__(self, key):
        super(OrderedDotDict, self).__delitem__(key)
        self.__keyorder__.__delitem__(key)

    def __delattr__(self, attr):
        super(OrderedDotDict, self).__delattr__(attr)
        self.__keyorder__.__delitem__(attr)

    def pop(self, *args, **kwargs):
        rv = super(OrderedDotDict, self).pop( *args, **kwargs)
        self.__keyorder__.pop(*args, **kwargs)
        return rv

    def update(self, *args, **kwargs):
        for k,v in collections.OrderedDict(*args, **kwargs).iteritems():
            self[k] = v


class SortedCollection(object):
    """Sequence sorted by a key function.

    TAKEN WITHOUT MODIFICATION FROM:
    https://code.activestate.com/recipes/577197-sortedcollection/
    (EDIT: removed ``__reduce__`` - better behavior with __dict__ states)

    SortedCollection() is much easier to work with than using bisect() directly.
    It supports key functions like those use in sorted(), min(), and max().
    The result of the key function call is saved so that keys can be searched
    efficiently.

    Instead of returning an insertion-point which can be hard to interpret, the
    five find-methods return a specific item in the sequence. They can scan for
    exact matches, the last item less-than-or-equal to a key, or the first item
    greater-than-or-equal to a key.

    Once found, an item's ordinal position can be located with the index() method.
    New items can be added with the insert() and insert_right() methods.
    Old items can be deleted with the remove() method.

    The usual sequence methods are provided to support indexing, slicing,
    length lookup, clearing, copying, forward and reverse iteration, contains
    checking, item counts, item removal, and a nice looking repr.

    Finding and indexing are O(log n) operations while iteration and insertion
    are O(n).  The initial sort is O(n log n).

    The key function is stored in the 'key' attibute for easy introspection or
    so that you can assign a new key function (triggering an automatic re-sort).

    In short, the class was designed to handle all of the common use cases for
    bisect but with a simpler API and support for key functions.

    >>> from pprint import pprint
    >>> from operator import itemgetter

    >>> s = SortedCollection(key=itemgetter(2))
    >>> for record in [
    ...         ('roger', 'young', 30),
    ...         ('angela', 'jones', 28),
    ...         ('bill', 'smith', 22),
    ...         ('david', 'thomas', 32)]:
    ...     s.insert(record)

    >>> pprint(list(s))         # show records sorted by age
    [('bill', 'smith', 22),
     ('angela', 'jones', 28),
     ('roger', 'young', 30),
     ('david', 'thomas', 32)]

    >>> s.find_le(29)           # find oldest person aged 29 or younger
    ('angela', 'jones', 28)
    >>> s.find_lt(28)           # find oldest person under 28
    ('bill', 'smith', 22)
    >>> s.find_gt(28)           # find youngest person over 28
    ('roger', 'young', 30)

    >>> r = s.find_ge(32)       # find youngest person aged 32 or older
    >>> s.index(r)              # get the index of their record
    3
    >>> s[3]                    # fetch the record at that index
    ('david', 'thomas', 32)

    >>> s.key = itemgetter(0)   # now sort by first name
    >>> pprint(list(s))
    [('angela', 'jones', 28),
     ('bill', 'smith', 22),
     ('david', 'thomas', 32),
     ('roger', 'young', 30)]
    """

    def __init__(self, iterable=(), key=None):
        self._given_key = key
        key = (lambda x: x) if key is None else key
        decorated = sorted((key(item), item) for item in iterable)
        self._keys = [k for k, item in decorated]
        self._items = [item for k, item in decorated]
        self._key = key

    def _getkey(self):
        return self._key

    def _setkey(self, key):
        if key is not self._key:
            self.__init__(self._items, key=key)

    def _delkey(self):
        self._setkey(None)

    key = property(_getkey, _setkey, _delkey, 'key function')

    def clear(self):
        self.__init__([], self._key)

    def copy(self):
        return self.__class__(self, self._key)

    def __len__(self):
        return len(self._items)

    def __getitem__(self, i):
        return self._items[i]

    def __iter__(self):
        return iter(self._items)

    def __reversed__(self):
        return reversed(self._items)

    def __repr__(self):
        return '%s(%r, key=%s)' % (
            self.__class__.__name__,
            self._items,
            getattr(self._given_key, '__name__', repr(self._given_key))
        )

    def __contains__(self, item):
        k = self._key(item)
        i = bisect_left(self._keys, k)
        j = bisect_right(self._keys, k)
        return item in self._items[i:j]

    def index(self, item):
        'Find the position of an item.  Raise ValueError if not found.'
        k = self._key(item)
        i = bisect_left(self._keys, k)
        j = bisect_right(self._keys, k)
        return self._items[i:j].index(item) + i

    def count(self, item):
        'Return number of occurrences of item'
        k = self._key(item)
        i = bisect_left(self._keys, k)
        j = bisect_right(self._keys, k)
        return self._items[i:j].count(item)

    def insert(self, item):
        'Insert a new item.  If equal keys are found, add to the left'
        k = self._key(item)
        i = bisect_left(self._keys, k)
        self._keys.insert(i, k)
        self._items.insert(i, item)

    def insert_right(self, item):
        'Insert a new item.  If equal keys are found, add to the right'
        k = self._key(item)
        i = bisect_right(self._keys, k)
        self._keys.insert(i, k)
        self._items.insert(i, item)

    def remove(self, item):
        'Remove first occurence of item.  Raise ValueError if not found'
        i = self.index(item)
        del self._keys[i]
        del self._items[i]

    def find(self, k):
        'Return first item with a key == k.  Raise ValueError if not found.'
        i = bisect_left(self._keys, k)
        if i != len(self) and self._keys[i] == k:
            return self._items[i]
        raise ValueError('No item found with key equal to: %r' % (k,))

    def find_le(self, k):
        'Return last item with a key <= k.  Raise ValueError if not found.'
        i = bisect_right(self._keys, k)
        if i:
            return self._items[i-1]
        raise ValueError('No item found with key at or below: %r' % (k,))

    def find_lt(self, k):
        'Return last item with a key < k.  Raise ValueError if not found.'
        i = bisect_left(self._keys, k)
        if i:
            return self._items[i-1]
        raise ValueError('No item found with key below: %r' % (k,))

    def find_ge(self, k):
        'Return first item with a key >= equal to k.  Raise ValueError if not found'
        i = bisect_left(self._keys, k)
        if i != len(self):
            return self._items[i]
        raise ValueError('No item found with key at or above: %r' % (k,))

    def find_gt(self, k):
        'Return first item with a key > k.  Raise ValueError if not found'
        i = bisect_right(self._keys, k)
        if i != len(self):
            return self._items[i]
        raise ValueError('No item found with key above: %r' % (k,))