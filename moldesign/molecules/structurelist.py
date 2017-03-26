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


import moldesign as mdt
from moldesign import utils
from . import AtomContainer, AtomList


class ChildList(AtomContainer):
    """ A list of biochemical objects that can be accessed by name or by index.
    """
    __len__ = utils.Alias('_childinorder.__len__')
    __iter__ = utils.Alias('_childinorder.__iter__')
    index = utils.Alias('_childinorder.index')

    def __str__(self):
        return str(self._childinorder._items)

    def __repr__(self):
        try:
            return '<Children of %s: %s>' % (self.parent, self)
        except (KeyError, AttributeError):
            return '<ChildList @ %x (exception in self.__repr__)>' % id(self)

    def __init__(self, parent, index_children=False):
        super(ChildList, self).__init__()
        self.parent = parent
        self._childbyname = {}
        self._childinorder = utils.SortedCollection(key=_sortkey)
        self.index_children = index_children

    def __dir__(self):
        return self.__dict__.keys() + self.__class__.__dict__.keys() + self._childbyname.keys()

    def __getitem__(self, item):
        if isinstance(item, int) and item < 0:
            item = len(self) + item
        if isinstance(item, basestring):
            if item not in self._childbyname:
                raise KeyError('No object in "%s" named "%s"' % (self.parent, item))
            return self._childbyname[item]
        else:
            try:
                return self._childinorder[item]
            except IndexError:
                raise IndexError("No object with index '%d' in %s" % (item, self.parent))

    def __setitem__(self, key, val):
        if key in self._childbyname:
            raise KeyError('%s already exists in %s' % (key, self.parent))
        self._childbyname[key] = val
        self._childinorder.insert_right(val)
        self.rebuild()

    def __contains__(self, item):
        if isinstance(item, basestring) or item is None:
            return (item in self._childbyname)
        else:
            return (item in self._childinorder)

    def _remove(self, item):
        self._childinorder.resort()  # this is a hack because the indices often change
        self._childinorder.remove(item)
        self.rebuild()

    def __getattr__(self, item):
        if item == '_childbyname':
            return self.__getattribute__('_childbyname')

        try:
            return self._childbyname[item]
        except KeyError:
            raise AttributeError('ChildList object in %s has no attribute %s.' %
                                 (self.parent.__repr__(), item))

    def iteratoms(self):
        """Iterate over all atoms

        Yields:
            Atom: all atoms in this entity and/or its children
        """
        for child in self:
            if isinstance(child, mdt.Atom):
                yield child
            else:
                for atom in child.iteratoms():
                    yield atom

    @property
    def atoms(self):
        """ AtomList: a sorted list of all atoms in this entity and/or its children
        """
        return AtomList(self.iteratoms())

    def rebuild(self):
        self._childbyname = {obj.name: obj for obj in self._childinorder}
        if self.index_children:
            for idx, child in enumerate(self._childinorder):
                child._index = idx
        self._childinorder.resort()



def _sortkey(x):
    return (x.index, x.pdbindex)