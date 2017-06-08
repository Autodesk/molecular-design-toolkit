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
from past.builtins import basestring
from sortedcontainers import SortedListWithKey
import functools

import moldesign as mdt
from .. import utils
from . import toplevel, AtomList, AtomContainer


class ChildList(AtomContainer):
    """ A list of biochemical objects that can be accessed by name or by index.
    """
    __len__ = utils.Alias('_childinorder.__len__')
    __iter__ = utils.Alias('_childinorder.__iter__')
    index = utils.Alias('_childinorder.index')

    def __str__(self):
        return str([item.name for item in self])

    def __repr__(self):
        try:
            return '<Children of %s: %s>' % (self.parent, self)
        except (KeyError, AttributeError):
            return '<ChildList @ %x (exception in self.__repr__)>' % id(self)

    def __init__(self, parent):
        super().__init__()
        self.parent = parent
        self._childbyname = {}
        self._childinorder = SortedListWithKey(key=_SortKey)

    def __dir__(self):
        return (list(self.__dict__.keys()) + list(self.__class__.__dict__.keys())
                + list(self._childbyname.keys()))

    def __getitem__(self, item):
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
        self._childinorder.add(val)

    def __contains__(self, item):
        if isinstance(item, basestring) or item is None:
            return (item in self._childbyname)
        else:
            try:
                return (item in self._childinorder)
            except TypeError:  # meaning that the comparison couldn't even be made
                return False

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
        self._childinorder = SortedListWithKey([obj for obj in self._childbyname.values()],
                                               key=_SortKey)


@functools.total_ordering
class _SortKey(object):
    def __init__(self, obj):
        self.obj = obj

    def __eq__(self, other):
        return self.obj is other.obj

    def __ne__(self, other):
        # if not present, this always returns True in python 2. (???)
        return not self == other

    def __lt__(self, other):
        try:
            return self.obj.index < other.obj.index
        except TypeError:
            pass

        try:
            return self.obj.pdbindex < other.obj.pdbindex
        except TypeError:
            return id(self.obj) < id(other.obj)


@toplevel
class BioContainer(AtomContainer):
    """
    Generalized storage mechanism for hierarchical representation of biomolecules,
    e.g. by residue, chain, etc. Permits other groupings, provided that everything is
    tree-like.

    All children of a given entity must have unique names. An individual child can be retrieved with
    ``biocontainer.childname`` or ``biocontainer['childname']`` or ``biocontainer[index]``

    Yields:
        BioContainer or mdt.Atom: this entity's children, in order
    """

    __getitem__ = utils.Alias('children.__getitem__')
    __len__ = utils.Alias('children.__len__')
    __iter__ = utils.Alias('children.__iter__')
    __contains__ = utils.Alias('children.__contains__')
    atoms = utils.Alias('children.atoms')
    iteratoms = utils.Alias('children.iteratoms')
    rebuild = utils.Alias('children.rebuild')

    def __init__(self, name=None, molecule=None, index=None, pdbname=None, pdbindex=None,
                 **kwargs):
        """  Initialization:

        Args:
            name (str): Name of this biocontainer
            parent (mdt.Molecule): molecule this biocontainer belongs to
            index (int): index of this biocontainer in the parent molecule
            pdbname (str): PDB-format name of this biocontainer
            pdbindex (str): Index of this biocontainer in PDB format
        """
        super().__init__()
        self.children = ChildList(self)
        self.molecule = molecule
        self.name = name
        self.index = index

        self.pdbname = pdbname
        self.pdbindex = pdbindex

        for name, val in kwargs.items():
            setattr(self, name, val)

    def add(self, item, key=None):
        """ Add a child to this entity.

        Raises:
            KeyError: if an object with this key already exists

        Args:
            item (BioContainer or mdt.Atom): the child object to add
            key (str): Key to retrieve this item (default: ``item.name`` )
        """
        if key is None:
            key = item.name
        self.children[key] = item

    __setitem__ = add

    def __dir__(self):
        return (list(self.__dict__.keys()) +
                list(self.__class__.__dict__.keys()) +
                [x.name for x in self.children])

    def __hash__(self):
        """ Explicitly hash by object id
        """
        return id(self)

    def __eq__(self, other):
        return self is other

    def __repr__(self):
        try:
            if self.molecule is not None:
                return '<%s in %s>' % (self, self.molecule)
            else:
                return '<%s (no molecule)>' % self
        except (KeyError, AttributeError):
            return '<%s at %s (exception in __repr__)>' % (self.__class__.__name__,
                                                     id(self))

    def __str__(self):
        return '%s %s (index=%s)' % (self.__class__.__name__,
                                     self.name, str(self.index))

    def __call__(self, **kwargs):
        """
        Allow for some simple queries, i.e. mol.chain['A'].residue(pdbname='ALA')
        """
        retlist = []
        for child in self:
            for key, val in kwargs.items():
                if hasattr(child, key) and getattr(child, key) == val:
                    retlist.append(child)
        return retlist


@toplevel
class Instance(BioContainer):
    """ The singleton biomolecular container for each ``Molecule``. Its children are generally
    PDB chains. Users won't ever really see this object.
    """
    def __str__(self):
        return str(self.children)

    def __repr__(self):
        return '<Molecule instance: %s>' % str(self.children)


