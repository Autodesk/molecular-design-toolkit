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

import moldesign as mdt
from .. import utils, data
from . import BioContainer, toplevel, HasResidues


@toplevel
class Chain(BioContainer, HasResidues):
    """ Biomolecular chain class - its children are almost always residues.

    Attributes:
        parent (mdt.Molecule): the molecule this residue belongs to
        chain (Chain): the chain this residue belongs to
    """
    @utils.args_from(BioContainer)
    def __init__(self, name=None, **kwargs):
        for key in ('pdbname', 'pdbindex'):
            val = kwargs.pop(key, None)
            if val is not None:
                if name is not None and val != name:
                    raise ValueError('Inconsistent name for chain: %s' % kwargs)
                name = val

        super().__init__(name=name, **kwargs)
        self._type = None

        self._5p_end = self._3p_end = self._n_terminal = self._c_terminal = None

    @property
    def pdbindex(self):
        return self.name

    @pdbindex.setter
    def pdbindex(self, val):
        if val is not None:
            self.name = val

    pdbname = pdbindex

    @property
    def type(self):
        """ str: the type of chain - protein, DNA, solvent, etc.

        This field returns the type of chain, classified by the following rules:
        1) If the chain contains only one type of residue, it is given that classification
           (so a chain containing only ions has type "ion"
        2) If the chain contains a biopolymer + ligands and solvent, it is classified as a
           biopolymer (i.e. 'protein', 'dna', or 'rna'). This is the most common case with .pdb
           files from the PDB.
        3) If the chain contains multiple biopolymer types, it will be given a hybrid classification
           (e.g. 'dna/rna', 'protein/dna') - this is rare!
        4) If it contains multiple kinds of non-biopolymer residues, it will be called "solvent"
           (if all non-bio residues are water/solvent/ion) or given a hybrid name as in 3)
        """
        if self._type is None:

            counts = collections.Counter(x.type for x in self.residues)
            unique_types = sum(bool(v) for v in counts.values())
            if unique_types == 1:
                if self.num_residues == 1:
                    self._type = data.CHAIN_MONOMER_NAMES.get(self.residues[0].type,
                                                              self.residues[0].type)
                else:
                    self._type = self.residues[0].type

            else:
                polymer_types = sum(bool(counts[t]) for t in data.BIOPOLYMER_TYPES)
                if polymer_types == 1:  # the most common case - a polymer + solvent/ligands
                    for residue in self.residues:
                        if residue.type in data.BIOPOLYMER_TYPES: break
                    else:
                        assert False, "No biopolymer found but polymer_types==1"
                    self._type = residue.type

                elif polymer_types > 1:  # for rare cases, e.g. "DNA/RNA/PROTEIN"
                    self._type = '/'.join(k for k in data.BIOPOLYMER_TYPES if counts[k])
                elif polymer_types == 0:
                    if counts['unknown'] > 0:  # some molecule + solvent
                        self._type = '/'.join(k for k in counts if counts[k])
                    else:  # just solvent
                        self._type = 'solvent'
        return self._type

    @property
    def polymer_residues(self):
        for res in self.residues:
            if res.type in ('dna', 'protein'):
                yield res

    @property
    def solvent_residues(self):
        for res in self.residues:
            if res.type in ('water', 'solvent', 'ion'):
                yield res

    @property
    def unclassified_residues(self):
        for res in self.residues:
            if res.type == 'unknown':
                yield res

    def get_ligand(self):
        """ Return a (single) ligand if it exists; raises ValueError if there's not exactly one

        This is a utility routine to get a single ligand from a chain. If there's exactly one
        residue, it is returned. If not, ValueError is raised - use
        :meth:`Chain.unclassified_residues` to get an iterator over all unclassified residues.

        Returns:
            moldesign.Residue: ligand residue

        Raises:
            ValueError: if the chain does not contain exactly one unclassifiable residue
        """
        iterator = self.unclassified_residues
        try:
            ligand = next(iterator)
        except StopIteration:
            raise ValueError('This chain does not appear to contain any ligands')

        try:
            nextligand = next(iterator)
        except StopIteration:
            return ligand
        else:
            raise ValueError('Multiple ligands detected. Use `chain.unclassified_residues` to '
                             'iterate over them')

    def copy(self):
        """ Create a copy of this chain and all topology within.

        Returns:
            Chain: copy of this chain and all atoms/residues it contains. This copy will NOT
            reference any parent molecule
        """
        newatoms = super().copy_atoms()
        return newatoms[0].chain

    @property
    def num_residues(self):
        return len(self)

    nresidues = numresidues = num_residues

    @property
    def residues(self):
        """ChildList: list of residues in this chain """
        return self.children

    def add(self, residue, **kwargs):
        if residue.chain is None:
            residue.chain = self
        else:
            assert residue.chain is self, "Residue is not a member of this chain"

        return super().add(residue, **kwargs)

    def _get_chain_end(self, restype, selfattr, test):
        currval = getattr(self, selfattr)
        if currval is None or not getattr(currval, test):
            for residue in self.residues:
                if residue.type != restype:
                    continue
                if getattr(residue, test):
                    setattr(self, selfattr, residue)
                    break
        return getattr(self, selfattr)

    @property
    def c_terminal(self):
        """moldesign.Residue: The chain's C-terminus (or ``None`` if it does not exist)"""
        return self._get_chain_end('protein', '_c_terminal', 'is_c_terminal')

    @property
    def n_terminal(self):
        """moldesign.Residue: The chain's N-terminus (or ``None`` if it does not exist)"""
        return self._get_chain_end('protein', '_n_terminal', 'is_n_terminal')

    @property
    def fiveprime_end(self):
        """moldesign.Residue: The chain's 5' base (or ``None`` if it does not exist)"""
        return self._get_chain_end('dna', '_5p_end', 'is_5prime_end')

    @property
    def threeprime_end(self):
        """moldesign.Residue: The chain's 3' base (or ``None`` if it does not exist)"""
        return self._get_chain_end('dna', '_3p_end', 'is_3prime_end')

    def assign_biopolymer_bonds(self):
        """Connect bonds between residues in this chain.

        See Also:
            :ref:`moldesign.Residue.assign_template_bonds`

        Raises:
            ValueError: if ``residue.resname`` is not in bioresidue templates
            KeyError: if an atom in this residue is not recognized """
        residues = list(self)
        residues.sort(key=lambda x: int(x.pdbindex))
        bond_graph = {}
        for ires in range(len(residues)-1):
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
        return self._get_sequence()

    def _get_sequence(self, _html=False):
        if _html:
            outputs = ['<span style="word-wrap:break-word;font-family:monospace">']
        else:
            outputs = []

        # get a list of ALL residues, including missing ones
        missing = self.molecule.metadata.get('missing_residues', {}).get(self.pdbname, {})
        all_residues = list(self.residues)
        all_residues.extend(mdt.helpers.MissingResidue(self.pdbname, resname, resnum)
                            for resnum, resname  in missing.items())
        all_residues.sort(key=lambda r: r.pdbindex)

        bracket_open = False
        for res in all_residues:
            if res.type not in ('protein', 'dna', 'rna'):  # don't display non-biological residues
                continue

            elif getattr(res, 'missing', False):  # deal with missing residues
                if bracket_open:
                    outputs.append(']')
                    bracket_open = False

                if _html:
                    outputs.append(
                        '<span title="%s%s (MISSING), chain %s" style="color:Red">%s</span>' %
                        (res.resname, res.pdbindex, self.pdbname, res.code))
                else:
                    outputs.append('-')

            elif res.code == '?':  # unknown residue type
                if not bracket_open:
                    outputs.append('[')
                    bracket_open = True
                else:
                    outputs.append(',')

                if _html:
                    outputs.append('<span title="%s (idx %d), chain %s">%s</span>' %
                                   (res.name, res.index, res.chain.name, res.pdbname))
                else:
                    outputs.append(res.pdbname)

            else:  # we have a normal residue name
                if bracket_open:
                    outputs.append(']')
                    bracket_open = False

                if _html:
                    outputs.append('<span title="%s (idx %d), chain %s">%s</span>' %
                                   (res.name, res.index, res.chain.name, res.code))
                else:
                    outputs.append(res.code)

        if bracket_open:
            outputs.append(']')

        if _html:
            outputs.append('</span>')

        return ''.join(outputs)

    def __str__(self):
        return 'Chain %s' % self.name
