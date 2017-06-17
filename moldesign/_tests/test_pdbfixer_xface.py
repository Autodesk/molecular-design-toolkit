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

import numpy as np

import moldesign as mdt
from moldesign import units as u

import pytest

from .molecule_fixtures import pdb3aid, benzene, pdb1yu8


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
    newmol = mdt.mutate_residues(mol, {mol.residues[residx]: resname})

    assert newmol.num_chains == mol.num_chains
    assert mol.num_residues == newmol.num_residues
    foundnew = False

    for i, (res, newres) in enumerate(zip(mol.residues, newmol.residues)):
        if i == residx:
            foundnew = True
            assert newres.resname == resname
            assert newres.name == resname+str(newres.pdbindex)
            atomnames = set(atom.name for atom in newres)
            assert len(atomnames) == newres.num_atoms
            assert atomnames.issubset(allatoms)
        else:
            assert res.name == newres.name
            assert res.num_atoms == newres.num_atoms

            for oldatom, newatom in zip(res, newres):
                assert oldatom.name == newatom.name
                assert oldatom.atnum == newatom.atnum
                if not foundnew:
                    assert oldatom.pdbindex == newatom.pdbindex


def test_mutate_docstring_dict_example(pdb3aid):
    mol = pdb3aid
    assert mol.residues[5].resname != 'ALA'
    mut = mdt.mutate_residues(mol, {mol.residues[5]: 'ALA'})  # mutate residue 5 to ALA
    assert mut.residues[5].resname == 'ALA'


def test_mutation_nomenclature_string_only(pdb3aid):
    mol = pdb3aid
    res25 = mol.get_residues(pdbindex=25)
    assert len(res25) == 2
    assert [r.resname for r in res25] == ['ASP', 'ASP']

    mut = mdt.mutate_residues(mol, 'D25M')  # Mutate ALA43 to MET43
    assert mut.get_residues()
    mut25 = mut.get_residues(pdbindex=25)
    assert len(mut25) == 2
    assert [r.resname for r in mut25] == ['MET', 'MET']


def test_multiple_mutations(pdb3aid):
    mol = pdb3aid
    mut = mdt.mutate_residues(mol, ['A.2S', 'B.3S'])  # Mutate Chain A res 2 and B 3 to SER
    assert [r.resname for r in mut.chains['A'].get_residues(pdbindex=2)] == ['SER']
    assert [r.resname for r in mut.chains['B'].get_residues(pdbindex=3)] == ['SER']

