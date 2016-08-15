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
import collections

import moldesign as mdt
from moldesign import data, utils

from . import toplevel, AtomList, AtomContainer


@toplevel
class Entity(AtomContainer, utils.DictLike):
    """
    Generalized storage mechanism for hierarchical representation of biomolecules,
    e.g. by residue, chain, etc. Permits other groupings, provided that everything is
    tree-like.

    All children of a given entity must have unique keys. An individual child can be retrieved with
    ``entity.childname`` or ``entity['childname']`` or ``entity[index]``

    Yields:
        Entity or mdt.Atom: this entity's children, in order
    """
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
        self.molecule = molecule
        self.name = name
        self.index = index

        self.pdbname = pdbname
        self.pdbindex = pdbindex

        self._indexlist = None

        for name, val in kwargs.iteritems():
            setattr(self, name, val)

    @property
    def mass(self):
        """ u.Scalar[mass]: total mass of this object
        """
        return sum(self.atoms.mass)

    def add(self, item, key=None):
        """ Add a child to this entity.

        Raises:
            KeyError: if an object with this key already exists

        Args:
            item (Entity or mdt.Atom): the child object to add
            key (str): Key to retrieve this item (default: ``item.name`` )
        """
        self._indexlist = None
        if key is None:
            key = item.name
        if key in self:
            raise KeyError('%s already exists in %s' % (key, self))
        self[key] = item

    def __contains__(self, item):
        # TODO: O(N) lookup is bad
        return item in self.children

    def __getitem__(self, item):
        try:
            return super(Entity, self).__getitem__(item)
        except (KeyError, TypeError) as orig:
            try:
                return self.indexlist[item]
            except (TypeError, IndexError) as exc:
                raise orig

    @property
    def indexlist(self):
        """ list: list of all children, in order
        """
        if self._indexlist is None:
            self._indexlist = list(self)
        return self._indexlist

    def __iter__(self):
        if self._indexlist is not None:
            for item in self._indexlist: yield item
        else:
            keys = sorted(self.keys())
            for k in keys: yield self[k]

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

    def iteratoms(self):
        """Iterate over all atoms (i.e. leaf nodes)

        Yields:
            Atom: all atoms in this entity and/or its children
        """
        for item in self.itervalues():
            if hasattr(item,'itervalues'):
                for x in item.itervalues(): yield x
            else:
                yield item

    @property
    def atoms(self):
        """ AtomList: a sorted list of all atoms in this entity and/or its children
        """
        atoms = AtomList(sorted(self.iteratoms(),
                                key=lambda x: x.index))
        return atoms



@toplevel
class Instance(Entity):
    """ The singleton biomolecular container for each ``Molecule``. Its children are generally
    PDB chains. Users won't ever really see this object.
    """
    def __str__(self):
        return "biounit container (chains: %s) for molecule %s" % (', '.join(self.keys()), self.molecule.name)

