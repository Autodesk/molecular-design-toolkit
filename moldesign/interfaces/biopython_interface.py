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

import itertools
import string

import Bio.PDB
import Bio.PDB.MMCIF2Dict
import numpy as np

import moldesign as mdt
from moldesign import units as u
from moldesign.helpers.pdb import BioAssembly
from moldesign.utils import exports


@exports
def biopython_to_mol(struc):
    """Convert a biopython PDB structure to an MDT molecule.

    Note:
        Biopython doesn't deal with bond data, so no bonds will be present
        in the Molecule

    Args:
        struc (Bio.PDB.Structure.Structure): Biopython PDB structure to convert

    Returns:
        moldesign.Molecule: converted molecule
    """
    # TODO: assign bonds using 1) CONECT records, 2) residue templates, 3) distance
    newatoms = []
    backup_chain_names = list(string.ascii_uppercase)

    for chain in struc.get_chains():
        tmp, pdbidx, pdbid = chain.get_full_id()
        if not pdbid.strip():
            pdbid = backup_chain_names.pop()
        newchain = mdt.Chain(pdbname=pdbid.strip())

        for residue in chain.get_residues():
            newresidue = mdt.Residue(pdbname=residue.resname.strip(),
                                     pdbindex=residue.id[1])

            newchain.add(newresidue)

            for atom in residue.get_atom():
                elem = atom.element
                if len(elem) == 2:
                    elem = elem[0] + elem[1].lower()
                newatom = mdt.Atom(element=elem,
                                   name=atom.get_name(),
                                   pdbname=atom.get_name(),
                                   pdbindex=atom.get_serial_number())
                newatom.position = atom.coord * u.angstrom
                newresidue.add(newatom)

                newatoms.append(newatom)

    return mdt.Molecule(newatoms, name=struc.get_full_id()[0])


def get_mmcif_assemblies(fileobj=None, mmcdata=None):
    """Parse an mmCIF file, return biomolecular assembly specifications

    Args:
        fileobj (file-like): File-like object for the PDB file
            (this object will be rewound before returning)
        mmcdata (dict): dict version of complete mmCIF data structure (if passed, this will
            not be read again from fileobj)

    Returns:
        Mapping[str, BioAssembly]: dict mapping assembly ids to BioAssembly instances
    """
    if mmcdata is None:
        mmcdata = get_mmcif_data(fileobj)

    if '_pdbx_struct_assembly.id' not in mmcdata:
        return {}  # no assemblies present

    # Get assembly metadata
    ids = mmcdata['_pdbx_struct_assembly.id']
    details = mmcdata['_pdbx_struct_assembly.details']
    chains = mmcdata['_pdbx_struct_assembly_gen.asym_id_list']
    opers = mmcdata['_pdbx_struct_assembly_gen.oper_expression']
    transform_ids = mmcdata['_pdbx_struct_oper_list.id']

    # Get matrix transformations
    tmat = np.zeros((4, 4)).tolist()
    for i in range(3):  # construct displacement vector
        tmat[i][3] = mmcdata['_pdbx_struct_oper_list.vector[%d]' % (i+1)]
    for i, j in itertools.product(range(0, 3), range(0, 3)):  # construct rotation matrix
        tmat[i][j] = mmcdata['_pdbx_struct_oper_list.matrix[%d][%d]' % (i+1, j+1)]
    transforms = _make_transform_dict(tmat, transform_ids)

    # Make sure it's a list
    if not isinstance(ids, list):
        ids = [ids]
        details = [details]
        chains = [chains]
        opers = [opers]

    # now create the assembly specifications
    assemblies = {}
    for id, detail, chainlist, operlist in zip(ids, details, chains, opers):
        assert id not in assemblies
        transforms = [transforms[i] for i in operlist.split(',')]
        assemblies[id] = BioAssembly(detail, chainlist.split(','), transforms)

    return assemblies


def _make_transform_dict(tmat, transform_ids):
    if isinstance(transform_ids, list):
        for i, j in itertools.product(range(0, 3), range(0, 4)):
            tmat[i][j] = list(map(float, tmat[i][j]))
        tmat[3][3] = [1.0]*len(transform_ids)
        tmat[3][0] = tmat[3][1] = tmat[3][2] = [0.0]*len(transform_ids)
        tmat = np.array(tmat)
        transforms = {id: tmat[:, :, i] for i, id in enumerate(transform_ids)}
    else:
        for i, j in itertools.product(range(0, 4), range(0, 4)):
            tmat[i][j] = float(tmat[i][j])
        tmat[3][3] = 1.0
        tmat = np.array(tmat)
        transforms = {transform_ids: tmat}

    return transforms


def get_mmcif_data(fileobj):
    mmcdata = Bio.PDB.MMCIF2Dict.MMCIF2Dict(fileobj)
    fileobj.seek(0)  # rewind for future access
    return mmcdata
