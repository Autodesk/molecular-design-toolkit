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
import operator

import moldesign as mdt
from moldesign import utils

from . import toplevel, AtomList, AtomContainer


class ChildList(AtomContainer):
    """ A list of biochemical objects that can be accessed by name or by index.
    """
    __len__ = utils.Alias('_childinorder.__len__')
    __iter__ = utils.Alias('_childinorder.__iter__')

    def __str__(self):
        return str(self._childinorder._items)

    def __repr__(self):
        try:
            return '<Children of %s: %s>' % (self.parent, self)
        except:
            return '<ChildList @ %x (__repr__ failed)>' % id(self)

    def __init__(self, parent):
        super(ChildList, self).__init__()
        self.parent = parent
        self._childbyname = {}
        self._childinorder = utils.SortedCollection(key=_sortkey)

    def __dir__(self):
        return self.__dict__.keys() + self.__class__.__dict__.keys() + self._childbyname.keys()

    def __getitem__(self, item):
        if isinstance(item, basestring):
            return self._childbyname[item]
        else:
            return self._childinorder[item]

    def __setitem__(self, key, val):
        if key in self._childbyname:
            raise KeyError('%s already exists in %s' % (key, self.parent))
        self._childbyname[key] = val
        self._childinorder.insert_right(val)

    def __contains__(self, item):
        if isinstance(item, basestring):
            return (item in self._childbyname)
        else:
            return (item in self._childinorder)

    def __getattr__(self, item):
        if not hasattr(self, '_childbyname'):
            raise AttributeError('Uninitialized')

        try:
            return self._childbyname[item]
        except KeyError:
            raise AttributeError('ChildList object in %s has no attribute %s.' % (
                self.parent, item))

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


def _sortkey(x):
    return x.pdbindex


@toplevel
class Entity(AtomContainer):
    """
    Generalized storage mechanism for hierarchical representation of biomolecules,
    e.g. by residue, chain, etc. Permits other groupings, provided that everything is
    tree-like.

    All children of a given entity must have unique names. An individual child can be retrieved with
    ``entity.childname`` or ``entity['childname']`` or ``entity[index]``

    Yields:
        Entity or mdt.Atom: this entity's children, in order
    """

    __getitem__ = utils.Alias('children.__getitem__')
    __len__ = utils.Alias('children.__len__')
    __iter__ = utils.Alias('children.__iter__')
    atoms = utils.Alias('children.atoms')
    iteratoms = utils.Alias('children.iteratoms')
    rebuild = utils.Alias('children.rebuild')

    def __init__(self, name=None, molecule=None, index=None, pdbname=None, pdbindex=None,
                 **kwargs):
        """  Initialization:

        Args:
            name (str): Name of this entity
            parent (mdt.Molecule): molecule this entity belongs to
            index (int): index of this entity in the parent molecule
            pdbname (str): PDB-format name of this entity
            pdbindex (str): Index of this entity in PDB format
        """
        super(Entity, self).__init__()
        self.children = ChildList(self)
        self.molecule = molecule
        self.name = name
        self.index = index

        self.pdbname = pdbname
        self.pdbindex = pdbindex

        for name, val in kwargs.iteritems():
            setattr(self, name, val)

    def add(self, item, key=None):
        """ Add a child to this entity.

        Raises:
            KeyError: if an object with this key already exists

        Args:
            item (Entity or mdt.Atom): the child object to add
            key (str): Key to retrieve this item (default: ``item.name`` )
        """
        if key is None:
            key = item.name
        self.children[key] = item

    __setitem__ = add

    def __dir__(self):
        return (self.__dict__.keys() +
                self.__class__.__dict__.keys() +
                [x.name for x in self])

    def __getattr__(self, item):
        if not hasattr(self, 'children'):
            raise AttributeError("Uninitialized")
        try:
            return self.children[item]
        except KeyError:
            raise AttributeError('%s has no attribute "%s"' % (self, item))

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
        except:
            return '<%s at %s (exc in __repr__)>' % (self.__class__.__name__,
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
            for key, val in kwargs.iteritems():
                if hasattr(child, key) and getattr(child, key) == val:
                    retlist.append(child)
        return retlist


@toplevel
class Instance(Entity):
    """ The singleton biomolecular container for each ``Molecule``. Its children are generally
    PDB chains. Users won't ever really see this object.
    """
    def __str__(self):
        return str(self.children)

    def __repr__(self):
        return '<Molecule instance: %s>' % str(self.children)


