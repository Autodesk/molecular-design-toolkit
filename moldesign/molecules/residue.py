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

import collections

import moldesign as mdt
from .. import utils, data
from . import BioContainer, AtomList, toplevel


@toplevel
class Residue(BioContainer):
    """ A biomolecular residue - most often an amino acid, a nucleic base, or a solvent
    molecule. In PDB structures, also often refers to non-biochemical molecules.

    Its children are almost always residues.

    Attributes:
        parent (mdt.Molecule): the molecule this residue belongs to
        chain (Chain): the chain this residue belongs to
    """

    def copy(self):
        """ Create a copy of this residue and all topology within.

        Returns:
            Residue: copy of this residue and all its atoms. This copy will NOT
            reference any parent molecule
        """
        newatoms = super().copy_atoms()
        return newatoms[0].residue

    @utils.args_from(BioContainer)
    def __init__(self, **kwargs):
        """ Initialization
        Args:
            **kwargs ():
        """
        self.chain = kwargs.get('chain', None)
        super().__init__(**kwargs)
        if self.index is None and self.molecule is not None:
            self.index = self.molecule.residues.index(self)
        self.chainindex = None
        self._backbone = None
        self._sidechain = None
        self._template_name = None
        self._name = None

    @property
    def name(self):
        if self._name is not None:
            return self._name
        elif self.pdbname is None:
            return None
        elif self.pdbindex is None:
            return self.pdbname
        elif self.pdbname[-1].isdigit():
            return '%s (seq # %s)' % (self.pdbname, self.pdbindex)
        else:
            return self.pdbname + str(self.pdbindex)

    def _repr_markdown_(self):
        return self.markdown_summary()

    def markdown_summary(self):
        """ Markdown-formatted information about this residue

        Returns:
            str: markdown-formatted string
        """
        if self.type == 'placeholder':
            return '`%s`'%repr(self)

        if self.molecule is None:
            lines = ["<h3>Residue %s</h3>"%self.name]
        else:
            lines = ["<h3>Residue %s (index %d)</h3>"%(self.name, self.index)]

        if self.type == 'protein':
            lines.append('**Residue codes**: %s / %s'%(self.resname, self.code))
        else:
            lines.append("**Residue code**: %s"%self.resname)
        lines.append('**Type**: %s'%self.type)
        if self.resname in data.RESIDUE_DESCRIPTIONS:
            lines.append('**Description**: %s'%data.RESIDUE_DESCRIPTIONS[self.resname])

        lines.append('**<p>Chain:** %s'%self.chain.name)

        lines.append('**PDB sequence #**: %d'%self.pdbindex)

        terminus = None
        if self.type == 'dna':
            if self.is_3prime_end:
                terminus = "3' end"
            elif self.is_5prime_end:
                terminus = "5' end"
        elif self.type == 'protein':
            if self.is_n_terminal:
                terminus = 'N-terminus'
            elif self.is_c_terminal:
                terminus = 'C-terminus'
        if terminus is not None:
            lines.append('**Terminal residue**: %s of chain %s'%(terminus, self.chain.name))

        if self.molecule is not None:
            lines.append("**Molecule**: %s"%self.molecule.name)

        lines.append("**<p>Number of atoms**: %s"%self.num_atoms)
        if self.backbone:
            lines.append("**Backbone atoms:** %s"%', '.join(x.name for x in self.backbone))
            lines.append("**Sidechain atoms:** %s"%', '.join(x.name for x in self.sidechain))
        else:
            lines.append("**Atom:** %s"%', '.join(x.name for x in self.atoms))

        return '<br>'.join(lines)

    @name.setter
    def name(self, value):
        self._name = value

    @property
    def atoms(self):
        return self.children

    @property
    def bonded_residues(self):
        """ List[moldesign.Residue]: list of residues this residue is bonded to
        """
        residues = collections.OrderedDict()
        for atom in self.bonded_atoms:
            residues[atom.residue] = None
        return list(residues.keys())

    def add(self, atom, key=None):
        """Deals with atom name clashes within a residue - common for small molecules"""
        if atom.residue is not None:
            assert atom.residue is self, 'Atom already assigned to a residue'
        atom.residue = self
        if atom.chain is None:
            atom.chain = self.chain
        else:
            assert atom.chain == self.chain, "Atom's chain does not match residue's chain"

        if key is not None or atom.name not in self.children:
            return super().add(atom, key=key)
        else:
            return super().add(atom, key='%s%s' % (atom.name, len(self)))
    add.__doc__ = BioContainer.add.__doc__

    @property
    def is_n_terminal(self):
        """bool: this is the last residue in a peptide

        Raises:
            ValueError: if this residue is not an amino acid
        """
        if self.type != 'protein':
            raise ValueError('%s is not a recognized peptide monomer' % self)
        return self._is_starting_residue

    @property
    def is_c_terminal(self):
        """bool: this is the first residue in a peptide

        Raises:
            ValueError: if this residue is not an amino acid
        """
        if self.type != 'protein':
            raise ValueError('%s is not a recognized peptide monomer' % self)
        return self._is_ending_residue

    @property
    def is_5prime_end(self):
        """bool: this is the first base in a strand

        Raises:
            ValueError: if this residue is not a DNA base
        """
        if self.type != 'dna':
            raise ValueError('%s is not a recognized nucleic acid monomer' % self)
        return self._is_starting_residue

    @property
    def is_3prime_end(self):
        """bool: this is the last base in a strand

        Raises:
            ValueError: if this residue is not a DNA base
        """
        if self.type != 'dna':
            raise ValueError('%s is not a recognized nucleic acid monomer' % self)
        return self._is_ending_residue

    @property
    def is_monomer(self):
        """bool: this residue is not part of a biopolymer
        """
        return self._is_ending_residue and self._is_starting_residue

    @property
    def _is_ending_residue(self):
        """bool: this is the last residue in a polymer"""
        try:
            nextres = self.next_residue
        except StopIteration:
            return True
        except KeyError:
            # If we're here, the residue is missing some atoms. We'll fall back to checking the
            # next residues in line
            if self.index == len(self.molecule.residues):
                return True
            else:
                print('WARNING: %s is missing expected atoms. Attempting to infer chain end' % \
                    self)
                nextres = self.molecule.residues[self.index + 1]
                return not self._same_polymer(nextres)

        else:
            return False

    @property
    def _is_starting_residue(self):
        """bool: this is the first residue in a polymer"""
        try:
            prevres = self.prev_residue
        except StopIteration:
            return True
        except KeyError:
            # If we're here, the residue is missing some atoms. We'll fall back to checking the
            # next residues in line
            if self.index <= 0:
                assert self.index == 0
                return True
            else:
                print('WARNING: %s is missing expected atoms. Attempting to infer chain start' % \
                    self)
                nextres = self.molecule.residues[self.index - 1]
                return not self._same_polymer(nextres)
        else:
            return False

    def _same_polymer(self, otherres):
        return (otherres.type == self.type) and (otherres.chain is self.chain)

    def assign_template_bonds(self):
        """Assign bonds from bioresidue templates.

        Only assigns bonds that are internal to this residue (does not connect different residues).
        The topologies here assume pH7.4 and may need to be corrected for other pHs

        See Also:
            :ref:`moldesign.Chain.assign_biopolymer_bonds` for assigning inter-residue bonds

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
            if len(self) == 1 and self.type not in ('water', 'ion'):
                print('INFO: no bonds assigned to residue %s' % self)
                return
            else:
                raise KeyError("No bonding template for residue '%s'" % resname)

        # intra-residue bonds
        bond_graph = {atom: {} for atom in self}
        for atom in self:
            for nbrname, order in bonds_by_name.get(atom.name, {}).items():
                try:
                    nbr = self[nbrname]
                except KeyError:  # missing atoms are normal (often hydrogen)
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
        if self.chain.type == 'protein':
            return self._get_neighbor('C', 'C-terminus')
        elif self.chain.type == 'dna':
            return self._get_neighbor("O3'", "3' end")
        else:
            raise NotImplementedError('We only deal with dna and amino acids right now')

    @property
    def prev_residue(self):
        """Residue: The next residue in the chain (in the N-direction for proteins, 5' direction for
            nucleic acids)

        Raises:
            NotImplementedError: If we don't know how to deal with this type of biopolymer
            StopIteration: If there isn't a previous residue (i.e. it's a 5'- or N-terminus)
        """
        if self.chain.type == 'protein':
            return self._get_neighbor('N', 'N-terminus')
        elif self.chain.type == 'dna':
            try:
                return self._get_neighbor("P", "5' end")
            except KeyError:
                raise StopIteration("No 5' bond for this residue")
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
        return data.RESIDUE_TYPES.get(self.resname, 'unknown')

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
                try:
                    self._backbone.append(self[name])
                except KeyError:
                    pass
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

    @property
    def is_standard_residue(self):
        """ bool: this residue is a "standard residue" for the purposes of a PDB entry.

        In PDB files, this will be stored using 'ATOM' if this is a standard residue
        and 'HETATM' records if not.

        Note:
            We currently define "standard" residues as those whose 3 letter residue code appears in
            the ``moldesign.data.RESIDUE_DESCRIPTIONS`` dictionary. Although this seems to work
            well, we'd welcome a PR with a less hacky method.

        References:
            PDB format guide: http://www.wwpdb.org/documentation/file-format
        """
        return self.resname in mdt.data.RESIDUE_DESCRIPTIONS

    def __str__(self):
        return 'Residue %s (index %s, chain %s)' % (self.name, self.index, self.chain.name)


class HasResidues(object):
    """ Mixin for classes that *contain* residues (i.e. Molecules and Chains)
    """

    def get_residues(self, **queries):
        """Allows keyword-based residue queries. Returns residues that match ALL queries.

        Args:
            **queries (dict): attributes (or residue attributes) to match

        Examples:
            >>> mol.get_residues(type='protein')  # returns all amino acid residues in molecule
            >>> mol.chains['A'].get_residues(resname='ALA')  # returns all alanines in chain A
            >>> mol.get_residues(chain='A')  # all residues in chain A


        Returns:
            List[Residues]: residues matching the query
        """
        if not queries:
            return list(self.residues)

        result = []
        for res in self.residues:
            for field, val in queries.items():
                if field == 'chain' and isinstance(val, basestring):
                    if res.chain.name != val:
                        break
                    else:
                        continue

                elif getattr(res, field, None) != val:
                    break
            else:
                result.append(res)

        return result