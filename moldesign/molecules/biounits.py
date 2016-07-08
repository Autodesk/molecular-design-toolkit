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

    @property
    def natoms(self):
        """ int: number of atoms contained in this entity and its children
        """
        return len(self.atoms)
    num_atoms = numatoms = natoms


@toplevel
class Instance(Entity):
    """ The singleton biomolecular container for each ``Molecule``. Its children are generally
    PDB chains. Users won't ever really see this object.
    """
    def __str__(self):
        return "biounit container (chains: %s) for molecule %s" % (', '.join(self.keys()), self.molecule.name)


@toplevel
class Residue(Entity):
    """ A biomolecular residue - most often an amino acid, a nucleic base, or a solvent
    molecule. In PDB structures, also often refers to non-biochemical molecules.

    Its children are almost always residues.

    Attributes:
        parent (mdt.Molecule): the molecule this residue belongs to
        chain (Chain): the chain this residue belongs to
    """

    def copy(self):
        newatoms = super(Residue, self).copy()
        return newatoms[0].residue
    copy.__doc__ = Entity.copy.__doc__

    @utils.args_from(Entity)
    def __init__(self, **kwargs):
        """ Initialization
        Args:
            **kwargs ():
        """
        self.chain = kwargs.get('chain', None)
        super(Residue, self).__init__(**kwargs)
        if self.index is None and self.molecule is not None:
            self.index = self.molecule.residues.index(self)
        self.chainindex = None
        self._backbone = None
        self._sidechain = None
        self._template_name = None
        if self.name is None: self.name = self.pdbname + str(self.pdbindex)

    def add(self, atom, key=None):
        """Deals with atom name clashes within a residue - common for small molecules"""
        if atom.residue is not None:
            assert atom.residue is self, 'Atom already assigned to a residue'
        atom.residue = self
        if atom.chain is None:
            atom.chain = self.chain
        else:
            assert atom.chain == self.chain, "Atom's chain does not match residue's chain"

        if key is not None or atom.name not in self:
            return super(Residue, self).add(atom, key=key)
        else:
            return super(Residue, self).add(atom, key='%s%s' % (atom.name, len(self)))
    add.__doc__ = Entity.__doc__

    @property
    def is_n_terminal(self):
        """bool: this is the last residue in a peptide

        Raises:
            ValueError: if this residue is not an amino acid
        """
        if self.type != 'protein':
            raise ValueError('This residue is not a recognized peptide monomer')
        return self._is_starting_residue

    @property
    def is_c_terminal(self):
        """bool: this is the first residue in a peptide

        Raises:
            ValueError: if this residue is not an amino acid
        """
        if self.type != 'protein':
            raise ValueError('This residue is not a recognized peptide monomer')
        return self._is_ending_residue

    @property
    def is_5prime_end(self):
        """bool: this is the first base in a strand

        Raises:
            ValueError: if this residue is not a DNA base
        """
        if self.type != 'dna':
            raise ValueError('This residue is not a recognized nucleic acid monomer')
        return self._is_starting_residue

    @property
    def is_3prime_end(self):
        """bool: this is the last base in a strand

        Raises:
            ValueError: if this residue is not a DNA base
        """
        if self.type != 'dna':
            raise ValueError('This residue is not a recognized nucleic acid monomer')
        return self._is_ending_residue

    @property
    def is_monomer(self):
        """bool: this residue is not part of a biopolymer
        """
        return self._is_ending_residue and self._is_starting_residue

    @property
    def _is_ending_residue(self):
        """bool: this is the last residue in the chain"""
        for otherres in self.molecule.residues[self.index+1:]:
            if otherres.chain != self.chain:
                return True
            elif otherres.type == self.type:
                return False
        return True

    @property
    def _is_starting_residue(self):
        """bool: this is the first residue of its type in the chain"""
        for otherres in self.molecule.residues[self.index-1::-1]:
            if otherres.chain != self.chain:
                return True
            elif otherres.type == self.type:
                return False

        return True

    def assign_template_bonds(self):
        """Assign bonds from bioresidue templates.

        Only assigns bonds that are internal to this residue (does not connect different residues).
        The topologies here assume pH7.4 and may need to be corrected for other pHs

        See Also:
            :ref:`moldesign.Chain.assign_biopolymer_bonds`

        Raises:
            ValueError: if ``residue.resname`` is not in bioresidue templates
            KeyError: if an atom in this residue is not recognized """
        try:
            resname = self.resname
            if self.type == 'protein':
                if self.is_n_terminal:
                    resname = self.resname + '_LSN3'  # the protonated form (with NH3+ on the end)
                elif self.is_c_terminal:
                    resname = self.resname + '_LEO2H'  # deprotonated form (COO-)
                elif self.is_monomer:
                    resname = self.resname + '_LFZW'  # free zwitterion form

            bonds_by_name = data.RESIDUE_BONDS[resname]
            self._template_name = resname
        except KeyError:
            raise ValueError("No bonding template for residue '%s'" % resname)

        # intra-residue bonds
        bond_graph = {atom: {} for atom in self}
        for atom in self:
            for nbrname, order in bonds_by_name.get(atom.name, {}).iteritems():
                try:
                    nbr = self[nbrname]
                except KeyError:  # missing atom in this structure is normal (often hydrogen)
                    pass
                else:
                    bond_graph[atom][nbr] = bond_graph[nbr][atom] = order

        # copy bonds into the right structure (do this last to avoid mangling the graph)
        for atom in bond_graph:
            atom.bond_graph.update(bond_graph[atom])

    @property
    def next_residue(self):
        """Residue:
            The next residue in the chain (in the C-direction for proteins, 3'
            direction for nucleic acids)

        Raises:
            NotImplementedError: If we don't know how to deal with this type of biopolymer
            StopIteration: If there isn't a next residue (i.e. it's a 3'- or C-terminus)
        """
        if self.type == 'protein':
            return self._get_neighbor('C', 'C-terminus')
        elif self.type == 'dna':
            return self._get_neighbor("O3'", "3' end")
        else:
            raise NotImplementedError('We only deal with dna and amino acids right now')

    @property
    def prev_residue(self):
        """Residue:
            The next residue in the chain (in the N-direction for proteins, 5' direction for nucleic acids)

        Raises:
            NotImplementedError: If we don't know how to deal with this type of biopolymer
            StopIteration: If there isn't a previous residue (i.e. it's a 5'- or N-terminus)
        """
        if self.type == 'protein':
            return self._get_neighbor('N', 'N-terminus')
        elif self.type == 'dna':
            return self._get_neighbor("P", "5' end")
        else:
            raise NotImplementedError('We only deal with dna and amino acids right now')

    def _get_neighbor(self, atomname, name):
        """Return the first residue found that's bound to the passed atom name
        """
        conn_atom = self[atomname]
        for nbr in conn_atom.bond_graph:
            if nbr.residue is not self:
                return nbr.residue
        else:
            raise StopIteration('%s reached' % name)

    @property
    def resname(self):
        """str: Synonym for pdbname"""
        return self.pdbname

    @resname.setter
    def resname(self, val):
        self.pdbname = val

    @property
    def type(self):
        """str: Classification of the residue (protein, solvent, dna, water, unknown)"""
        for restype, names in data.RESTYPES.iteritems():
            if self.pdbname in names:
                return restype
        return 'unknown'

    @property
    def code(self):
        """str: one-letter amino acid code or two letter nucleic acid code, or '?' otherwise"""
        return data.RESIDUE_ONE_LETTER.get(self.pdbname, '?')

    @property
    def atomnames(self):
        """Residue: synonym for ```self``` for for the sake of readability:
            ```molecule.chains['A'].residues[123].atomnames['CA']```
        """
        return self

    @property
    def backbone(self):
        """ AtomList: all backbone atoms for nucleic and protein residues
            (indentified using PDB names); returns None for other residue types
        """
        if self._backbone is None:
            if self.type not in data.BACKBONES:
                return None
            self._backbone = AtomList()
            for name in data.BACKBONES[self.type]:
                try: self._backbone.append(self[name])
                except KeyError: pass
        return self._backbone

    @property
    def sidechain(self):
        """ AtomList: all sidechain atoms for nucleic and protein residues
            (defined as non-backbone atoms); returns None for other residue types
        """
        if self._sidechain is None:
            if self.backbone is None:
                return None
            bb = set(self.backbone)
            self._sidechain = [atom for atom in self if atom not in bb]
        return self._sidechain

    def __str__(self):
        return 'Residue %s (index %d, chain %s)' % (self.name, self.index,
                                                    self.chain.name)

    def _repr_markdown_(self):
        return self.markdown_summary()

    def markdown_summary(self):
        """ Markdown-formatted information about this residue

        Returns:
            str: markdown-formatted string
        """
        if self.molecule is None:
            lines = ["<h3>Residue %s</h3>" % self.name]
        else:
            lines = ["<h3>Residue %s (index %d)</h3>" % (self.name, self.index)]

        if self.type == 'protein':
            lines.append('**Residue codes**: %s / %s' % (self.resname, self.code))
        else:
            lines.append("**Residue code**: %s" % self.resname)
        lines.append('**Type**: %s' % self.type)
        if self.resname in data.RESIDUE_DESCRIPTIONS:
            lines.append('**Description**: %s'%data.RESIDUE_DESCRIPTIONS[self.resname])

        lines.append('**<p>Chain:** %s' % self.chain.name)
        lines.append('**Sequence number**: %d' % self.pdbindex)
        if self.molecule is not None:
            lines.append("**Molecule**: %s" % self.molecule.name)

        lines.append("**<p>Number of atoms**: %s" % self.num_atoms)
        if self.backbone:
            lines.append("**Backbone atoms:** %s" % ', '.join(x.name for x in self.backbone))
            lines.append("**Sidechain atoms:** %s" % ', '.join(x.name for x in self.sidechain))
        else:
            lines.append("**Atom:** %s" % ', '.join(x.name for x in self.atoms))

        return '<br>'.join(lines)


@toplevel
class Chain(Entity):
    """ Biomolecular chain class - its children are almost always residues.

    Attributes:
        parent (mdt.Molecule): the molecule this residue belongs to
        chain (Chain): the chain this residue belongs to
    """
    @utils.args_from(Entity)
    def __init__(self, pdbname=None, **kwargs):
        super(Chain, self).__init__(pdbname=pdbname, **kwargs)
        if self.name is None: self.name = self.pdbname
        if self.pdbindex is not None: self.pdbindex = self.pdbname

    def copy(self):
        newatoms = super(Chain, self).copy()
        return newatoms[0].chain
    copy.__doc__ = Entity.copy.__doc__

    @property
    def residues(self):
        """Chain: synonym for 'self' to enhance readability,
            e.g. ```molecule.chains['A'].residue[123]```"""
        return self

    def add(self, residue, **kwargs):
        if residue.chain is None:
            residue.chain = self
        else:
            assert residue.chain is self, "Residue is not a member of this chain"

        return super(Chain, self).add(residue, **kwargs)

    def assign_biopolymer_bonds(self):
        """Connect bonds between residues in this chain. Does not .

        See Also:
            :ref:`moldesign.Residue.assign_template_bonds`

        Raises:
            ValueError: if ``residue.resname`` is not in bioresidue templates
            KeyError: if an atom in this residue is not recognized """
        residues = list(self)
        residues.sort(key=lambda x: int(x.pdbindex))
        bond_graph = {}
        for ires in xrange(len(residues)-1):
            r1 = residues[ires]
            r2 = residues[ires+1]

            # don't assign bonds unless these are contiguous bioresidues
            if r1.pdbindex + 1 != r2.pdbindex:
                continue
            restype = r1.type
            if r2.type != restype:
                continue

            # Create the bonds
            if restype == 'protein':
                bond_graph[r1['C']] = {r2['N']: 1}
                bond_graph[r2['N']] = {r1['C']: 1}
            elif restype == 'dna':
                bond_graph[r1["O3'"]] = {r2['P']: 1}
                bond_graph[r2['P']] = {r1["O3'"]: 1}
            elif restype == 'rna':
                raise NotImplementedError('RNA not yet implemented')

        # copy bonds into the right structure (do this last to avoid mangling the graph)
        for atom in bond_graph:
            atom.bond_graph.update(bond_graph[atom])

    @property
    def sequence(self):
        """str: this chain's residue sequence with one-letter residue codes
        """
        missing = '.'  # don't do this
        outputs = []
        last_idx = None
        for res in sorted(self, key=lambda x:x.pdbindex):
            if res.type not in ('protein', 'dna', 'rna'): continue
            if last_idx is not None:
                num_missing = res.pdbindex - last_idx - 1
                if num_missing > 0:
                    outputs.append(missing * (res.pdbindex - last_idx - 1))
            if res.code != '?':
                outputs.append(res.code)
            else:
                if len(outputs) > 0 and outputs[-1][-1] != ',':
                    outputs.append(',')
                outputs.append(res.pdbname + ',')
            last_idx = res.pdbindex
        return ''.join(outputs)

    def __str__(self):
        return 'Chain %s' % self.name
