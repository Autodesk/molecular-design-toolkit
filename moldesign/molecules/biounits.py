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

from . import toplevel, AtomContainer, ChildList


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

    INDEX_CHILDREN = False

    __getitem__ = utils.Alias('children.__getitem__')
    __len__ = utils.Alias('children.__len__')
    __iter__ = utils.Alias('children.__iter__')
    __contains__ = utils.Alias('children.__contains__')
    atoms = utils.Alias('children.atoms')
    iteratoms = utils.Alias('children.iteratoms')
    _rebuild = utils.Alias('children._rebuild')
    _remove = utils.Alias('children._remove')
    _renamechild = utils.Alias('children._rename')

    def __init__(self, name):
        """  Initialization:

        Args:
            name (str): Name of this biocontainer
            parent (mdt.Molecule): molecule this biocontainer belongs to
            index (int): index of this biocontainer in the parent molecule
            pdbname (str): PDB-format name of this biocontainer
            pdbindex (str): Index of this biocontainer in PDB format
        """
        super(BioContainer, self).__init__()
        self._index = None
        self._molecule = None
        self.name = name
        self.children = ChildList(self, self.INDEX_CHILDREN)

    def add(self, item):
        """ Add a child to this entity.

        Raises:
            KeyError: if an object with this key already exists

        Args:
            item (BioContainer or mdt.Atom): the child object to add
            key (str): Key to retrieve this item (default: ``item.name`` )
        """
        if item.name in self:
            item._name = self._getuniquename(item.name)
        assert item._name == item.name
        assert item.name not in self.children
        self.children[item.name] = item

    __setitem__ = add

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, n):
        self._name = n

    def _getuniquename(self, key):
        if key is None:
            key = ''
        assert isinstance(key, basestring)
        ix = 0
        k = key
        while k in self:
            ix += 1
            k = "%s.%d" % (key, ix)
        return k

    @property
    def index(self):
        return self._index

    @property
    def molecule(self):
        return self._molecule

    @molecule.setter
    def molecule(self, mol):
        if self.molecule is mol:
                return
        elif self.molecule is not None:
            raise ValueError("This object '%s' is already owned by %s" %
                             (self, self.molecule))
        else:
            self._molecule = mol
            for atom in self.atoms:
                atom.molecule = mol

    def __dir__(self):
        return (self.__dict__.keys() +
                self.__class__.__dict__.keys() +
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
            for key, val in kwargs.iteritems():
                if hasattr(child, key) and getattr(child, key) == val:
                    retlist.append(child)
        return retlist


@toplevel
class Instance(BioContainer):
    """ The singleton primary structure container for each ``Molecule``. Its children are generally
    PDB chains. Users won't ever really see this object.
    """
    INDEX_CHILDREN = True

    def __init__(self, molecule):
        """  Initialization:

        Args:
            molecule (moldesign.Molecule): the molecule this Instance represents
            name (str): Name of this instance
        """
        super(Instance, self).__init__(None)
        self._molecule = molecule

    @property
    def molecule(self):
        return self._molecule

    @property
    def name(self):
        return self.molecule.name + 'topology'

    @name.setter
    def name(self, value):
        assert value is None  # you don't actually set this name

    def __str__(self):
        return str(self.children)

    def __repr__(self):
        return '<Molecule instance: %s>' % str(self.children)

    def add(self, chain, _addatoms=True):
        utils.AutoIndexList.extend(self.molecule.residues, list(chain.residues))
        if _addatoms:
            utils.AutoIndexList.extend(self.molecule.atoms, list(chain.atoms))

        # rebuild all indices after adding children to a molecule
        self.children.rebuild()
        for c in self.children:
            c.children.rebuild()
            for residue in c.children:
                residue.children.rebuild()

        super(Instance, self).add(chain)
        chain._molecule = self.molecule

