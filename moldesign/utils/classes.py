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
        return self._keys.itervalues()

    def __len__(self):
        return len(self._keys)

    def __getitem__(self, item):
        return self._keys.values()[item]

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
            k = self._keys.keys()[index]
            return self._keys.pop(k)

    def __repr__(self):
        return '%s(%s)' % (type(self).__name__, self._keys.values())

    __str__ = __repr__


class DotDict(dict):
    """Dict with items accessible as attributes"""
    def __getstate__(self):
        retval = dict(__dict__=self.__dict__.copy(),
                      items=self.items())
        return retval

    def __setstate__(self, d):
        self.__dict__.update(d['__dict__'])
        self.update(d['items'])

    def __getattr__(self, item):
        try:
            return self[item]
        except KeyError:
            raise AttributeError("'%s' object has no attribute '%s'"
                                 % (self.__class__.__name__, item))

    def __setattr__(self, item, val):
        self[item] = val

    def __dir__(self):
        return dir(self.__class__) + self.keys()


class OrderedDotDict(collections.OrderedDict):
    """Dict with items accessible as attributes"""
    def __getstate__(self):
        retval = dict(__dict__=self.__dict__.copy(),
                      items=self.items())
        return retval

    def __setstate__(self, d):
        self.__dict__.update(d['__dict__'])
        self.update(d['items'])

    def __getattr__(self, item):
        try:
            return self[item]
        except KeyError:
            raise AttributeError("'%s' object has no attribute '%s'"
                                 % (self.__class__.__name__, item))

    def __setattr__(self, item, val):
        if item.startswith('_OrderedDict__') or item.startswith('__'):
            self.__dict__[item] = val
        else:
            self[item] = val

    def __dir__(self):
        return dir(self.__class__) + self.keys()


class Alias(object):
    """
    Descriptor that calls a child object's method.
    e.g.
    >>> class A(object):
    >>>     childkeys = Alias('child.keys')
    >>>     child = dict()
    >>>
    >>> a = A()
    >>> a.child['key'] = 'value'
    >>> a.childkeys() #calls a.child.keys(), returns ['key']
    ['key']
    """
    def __init__(self, objmethod):
        objname, methodname = objmethod.split('.')
        self.objname = objname
        self.methodname = methodname

    def __get__(self, instance, owner):
        proxied = getattr(instance, self.objname)
        return getattr(proxied,self.methodname)


class Synonym(object):
    """ An attribute (class or intance) that is just a synonym for another.
    """
    def __init__(self, name):
        self.name = name

    def __get__(self, instance, owner):
        return getattr(instance, self.name)

    def __set__(self, instance, value):
        setattr(instance, self.name, value)


class DictLike(object):
    """
    This just wraps normal dicts so that other classes don't have to inherit from a built-in class,
    which apparently breaks pickle quite frequently.
    """
    def __getstate__(self):
        retval = dict(__dict__=self.__dict__.copy())
        return retval

    def __setstate__(self, d):
        self.__dict__.update(d['__dict__'])

    def __init__(self, **kwargs):
        self.children = {}
        self.children.update(kwargs)

    def __getattr__(self, k):
        return getattr(self.children, k)

    __setitem__ = Alias('children.__setitem__')
    __getitem__ = Alias('children.__getitem__')

    def __len__(self):
        return len(self.children)


class Attribute(object):
    """For overriding a property in a superclass - turns the attribute back
    into a normal instance attribute"""
    def __init__(self, name):
        self.name = name

    def __get__(self, instance, cls):
        return getattr(instance, self.name)

    def __set__(self, instance, value):
        return setattr(instance, self.name, value)


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