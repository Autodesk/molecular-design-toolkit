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

from . import MolecularHierarchy, Residue, Chain, ChildList
from .. import utils

# TODO: fully delegate all biochemical property management and views to this object.
# TODO: atom.residue, atom.chain, mol.residues should all delegate to PrimaryStructure


class PrimaryStructure(MolecularHierarchy):
    """ A tree-like structure that encapsulates the biochemical organization of a ``Molecule``

    *Most* of the state related to the biochemical hierarchy is stored
    here, with one major exception: each atom stores a reference to its containing residue. That
    data must be consistent with the data in this object

    Args:
        mol (moldesign.Molecule): molecule this object stores hierarchy for
    """
    def __init__(self, mol):
        super().__init__(molecule=mol, name='root')
        self._defchain = Chain(name=None,
                               index=None,
                               molecule=None)
        self._defres = Residue(name=None,
                               index=None,
                               pdbindex=None,
                               pdbname=None,
                               chain=self._defchain,
                               molecule=None)
        self.residues = None

    def __str__(self):
        return str(self.children)

    def __repr__(self):
        return '<Molecule instance: %s>' % str(self.children)

    def rebuild_hierarchy(self):
        """
        Set up the chain/residue/atom hierarchy

        This rebuilds the primary structure tree using each object's parents. If atoms
        don't have a residue, they are assigned a default residue. If residues don't have a
        chain, they're assigned a default chain.
        """
        self.children = ChildList(self)

        conflicts = set()
        self._fix_atom_pdbnames(conflicts)
        self.residues = self._rebuild_residues()
        self._rebuild_chains(conflicts)

        if conflicts:
            print('WARNING: %s modified due to name clashes' % (', '.join(conflicts)))

    def _fix_atom_pdbnames(self, conflicts):
        """ Rename atom pdbindices to avoid name conflicts
        """
        pdbatoms = [atom for atom in self.molecule.atoms if atom.pdbindex is not None]
        if pdbatoms:
            min_pdb_atom = min(pdbatoms, key=lambda x: x.pdbindex)
            last_pdb_idx = min_pdb_atom.pdbindex-min_pdb_atom.index-1
        else:
            last_pdb_idx = 0

        for atom in self.molecule.atoms:
            # assign PDB indices
            if atom.pdbindex is None:
                atom.pdbindex = last_pdb_idx+1
            if last_pdb_idx is not None and atom.pdbindex <= last_pdb_idx:
                atom.pdbindex = last_pdb_idx+1
                conflicts.add('atom indices')
            last_pdb_idx = atom.pdbindex

    def _rebuild_residues(self):
        """ Recreate list of residues
        """
        num_biores = 0

        foundresidues = collections.OrderedDict()  # used as an ordered set
        for atom in self.molecule.atoms:
            res = atom.residue if atom.residue is not None else self._defres

            if res not in foundresidues:
                # rebuild first - names/ordering may have changed
                res.rebuild()

                # add it to the structure
                res.index = len(foundresidues)
                foundresidues[res] = None
                if res.type in ('dna', 'rna', 'protein'):
                    num_biores += 1

                if res.molecule is None:
                    res.molecule = self.molecule
                else:
                    assert res.molecule is self.molecule

            if atom not in res:
                res.add(atom)

        self.molecule.is_biomolecule = (num_biores >= 2)
        return list(foundresidues)

    def _rebuild_chains(self, conflicts):
        """ Recreate list of chains, renaming as necessary to avoid conflicts
        """
        for residue in self.residues:
            chain = residue.chain if residue.chain is not None else self._defchain

            if chain not in self:
                # rebuild first - names/ordering may have changed
                chain.rebuild()

                # rename chain to avoid conflicts if necessary
                oldname = chain.name
                if chain.name is None and len(self.residues) > 1:
                    chain.name = 'A'
                while chain.name in self:
                    if chain.name is None:
                        chain.name = 'A'
                    chain.name = chr(ord(chain.name)+1)

                # add it to the structure
                chain.index = len(self)
                self.add(chain)
                if chain.name != oldname:
                    conflicts.add('chain ids')

                if chain.molecule is None:
                    chain.molecule = self.molecule
                else:
                    assert chain.molecule is self.molecule

            if residue not in chain:
                chain.add(residue)

    def _repr_markdown_(self):
        """A markdown description of biomolecular structure.

        Returns:
            str: Markdown string"""
        lines = []

        if len(self.residues) > 1:
            table = self._get_residue_table()
            lines.append('### Residues')
            lines.append(table.markdown(replace={0: ' '})+'|')  # extra '|' is bug workaround (?)

            lines.append('### Biopolymer chains')
            seqs = []
            for chain in self:
                seq = chain._get_sequence(_html=True)
                if not seq.strip():  # don't write anything if there's no sequence
                    continue

                seqs.append('**%s**: %s' % (chain.name, seq))
            lines.append('<br>'.join(seqs))

        return lines

    def _get_residue_table(self):
        """Creates a data table summarizing this molecule's primary structure.

        Returns:
            moldesign.utils.MarkdownTable: markdown-formatted table"""
        table = utils.MarkdownTable(*(['chain']+
                                      'protein dna rna unknown water solvent ion'.split()))
        for chain in self:
            counts = {}
            unk = []
            for residue in chain.residues:
                cat = residue.type
                if cat == 'unknown':
                    unk.append(residue.name)
                counts[cat] = counts.get(cat, 0)+1
            counts['chain'] = '<b>%s</b>'%chain.name
            if 0 < len(unk) <= 4:
                counts['unknown'] = ','.join(unk)
            table.add_line(counts)
        return table

