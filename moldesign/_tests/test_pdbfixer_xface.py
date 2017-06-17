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

import moldesign as mdt
from .molecule_fixtures import pdb3aid


def test_inchain_residue_mutation_in_protein(pdb3aid):
    _mutate_and_check(pdb3aid, 3, 'ALA',
                      {'C', 'CA', 'CB', 'HA', 'HB1', 'HB2', 'HB3', 'HN2', 'N', 'O'})


def test_nterm_residue_mutate_protein(pdb3aid):
    _mutate_and_check(pdb3aid, 0, 'MET',
                      {'CE', 'HE1', 'C', 'N', 'HA', 'HB2', 'HE3', 'HG2', 'CG', 'CA', 'HG3',
                       'O', 'HB3', 'HN2', 'CB', 'HE2'})


def test_cterm_residue_mutate_protein(pdb3aid):
    cterm = pdb3aid.chains['A'].c_terminal
    _mutate_and_check(pdb3aid, cterm.index, 'LEU',
                      {'HD11', 'HD23', 'C', 'N', 'HD22', 'HA', 'HB2', 'CG', 'CA', 'CD2', 'HD21',
                       'O', 'CD1', 'HD12', 'HB3', 'HN2', 'HD13', 'CB', 'HG'})


def _mutate_and_check(mol, residx, resname, allatoms):
    newmol = mdt.mutate(mol, {mol.residues[residx]: resname})

    assert newmol.num_chains == mol.num_chains
    assert mol.num_residues == newmol.num_residues

    for i, (res, newres) in enumerate(zip(mol.residues, newmol.residues)):
        if i == residx:
            assert newres.resname == resname
            assert newres.name == resname+str(newres.pdbindex)
            atomnames = set(atom.name for atom in newres)
            assert len(atomnames) == newres.num_atoms
            assert atomnames.issubset(allatoms)
        else:
            assert res.name == newres.name
            assert res.num_atoms == newres.num_atoms

