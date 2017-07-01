from __future__ import print_function, absolute_import, division

import io

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

import itertools
import parmed
import future.utils

import moldesign as mdt
from .. import units as u
from .. import utils

if future.utils.PY2:
    from cStringIO import StringIO
else:
    from io import StringIO

def read_mmcif(f, reassign_chains=True):
    """Parse an mmCIF file (using the parmEd parser) and return a molecule

    Args:
        f (file): file-like object containing the mmCIF file
        reassign_chains (bool): reassign chain IDs from ``auth_asym_id`` to ``label_asym_id``

    Returns:
        moldesign.Molecule: parsed molecule
    """
    parmedmol = parmed.read_CIF(f)
    mol = parmed_to_mdt(parmedmol)
    if reassign_chains:
        f.seek(0)
        mol = _reassign_chains(f, mol)
    mdt.helpers.assign_biopolymer_bonds(mol)
    return mol


def read_pdb(f):
    """Parse an mmCIF file (using the parmEd parser) and return a molecule

    Args:
        f (file): file-like object containing the mmCIF file

    Returns:
        moldesign.Molecule: parsed molecule
    """
    parmedmol = parmed.read_PDB(f)
    mol = parmed_to_mdt(parmedmol)
    return mol


def write_pdb(mol, fileobj):
    """ Write a PDB file to a buffer

    Args:
        mol (moldesign.Molecule): molecule to write as pdb
        fileobj (io.IOBase): buffer to write to - bytes and text interfaces are acceptable
    """
    pmedmol = mol_to_parmed(mol)

    tempfile = StringIO()
    pmedmol.write_pdb(tempfile, renumber=False)

    if not isinstance(fileobj, io.TextIOBase) or 'b' in getattr(fileobj, 'mode', ''):
        binaryobj = fileobj
        fileobj = io.TextIOWrapper(binaryobj)
        wrapped = True
    else:
        wrapped = False

    _insert_conect_records(mol, pmedmol, tempfile, write_to=fileobj)

    if wrapped:
        fileobj.flush()
        fileobj.detach()

CONECT = 'CONECT %4d'

def _insert_conect_records(mol, pmdmol, pdbfile, write_to=None):
    """ Inserts TER records to indicate the end of the biopolymeric part of a chain

    Put CONECT records at the end of a pdb file that doesn't have them

    Args:
        mol (moldesign.Molecule): the MDT version of the molecule that pdbfile describes
        pdbfile (TextIO OR str): pdb file (file-like or string)

    Returns:
        TextIO OR str: copy of the pdb file with added TER records - it will be
         returned as the same type passed (i.e., as a filelike buffer or as a string)
    """
    conect_bonds = mdt.helpers.get_conect_pairs(mol)

    def get_atomseq(atom):
        return pmdmol.atoms[atom.index].number

    pairs = collections.OrderedDict()
    for atom, nbrs in conect_bonds.items():
        pairs[get_atomseq(atom)] = list(map(get_atomseq, nbrs))

    if isinstance(pdbfile, basestring):
        pdbfile = StringIO(pdbfile)

    if write_to is None:
        newf = StringIO()
    else:
        newf = write_to
    pdbfile.seek(0)

    for line in pdbfile:
        if line.split() == ['END']:
            for a1idx in pairs:
                for istart in range(0, len(pairs[a1idx]), 4):
                    pairindices = ''.join(("%5d" % idx) for idx in pairs[a1idx][istart:istart+4])
                    newf.write(str(CONECT % a1idx + pairindices + '\n'))

        newf.write(str(line))



def write_mmcif(mol, fileobj):
    mol_to_parmed(mol).write_cif(fileobj)


@utils.exports
def parmed_to_mdt(pmdmol):
    """ Convert parmed Structure to MDT Structure

    Args:
        pmdmol (parmed.Structure): parmed structure to convert

    Returns:
        mdt.Molecule: converted molecule
    """
    atoms = collections.OrderedDict()
    residues = {}
    chains = {}

    masses = [pa.mass for pa in pmdmol.atoms] * u.dalton
    positions = [[pa.xx, pa.xy, pa.xz] for pa in pmdmol.atoms] * u.angstrom

    for iatom, patm in enumerate(pmdmol.atoms):
        if patm.residue.chain not in chains:
            chains[patm.residue.chain] = mdt.Chain(pdbname=patm.residue.chain)
        chain = chains[patm.residue.chain]

        if patm.residue not in residues:
            residues[patm.residue] = mdt.Residue(resname=patm.residue.name,
                                                 pdbindex=patm.residue.number)
            residues[patm.residue].chain = chain
            chain.add(residues[patm.residue])
        residue = residues[patm.residue]

        atom = mdt.Atom(name=patm.name,
                        atnum=patm.atomic_number,
                        pdbindex=patm.number,
                        mass=masses[iatom])
        atom.position = positions[iatom]

        atom.residue = residue
        residue.add(atom)
        assert patm not in atoms
        atoms[patm] = atom

    for pbnd in pmdmol.bonds:
        atoms[pbnd.atom1].bond_to(atoms[pbnd.atom2], int(pbnd.order))

    mol = mdt.Molecule(list(atoms.values()),
                       metadata=_get_pdb_metadata(pmdmol))
    return mol


def _get_pdb_metadata(pmdmol):
    metadata = utils.DotDict(description=pmdmol.title)

    authors = getattr(pmdmol, 'journal_authors', None)
    if authors:
        metadata.pdb_authors = authors

    experimental = getattr(pmdmol, 'experimental', None)
    if experimental:
        metadata.pdb_experimental = experimental

    box_vectors = getattr(pmdmol, 'box_vectors', None)
    if box_vectors:
        metadata.pdb_box_vectors = box_vectors

    doi = getattr(pmdmol, 'doi', None)
    if doi:
        metadata.pdb_doi = doi
        metadata.url = "http://dx.doi.org/%s" % doi

    return metadata


@utils.exports
def mol_to_parmed(mol):
    """ Convert MDT Molecule to parmed Structure
    Args:
        mol (moldesign.Molecule):

    Returns:
        parmed.Structure
    """
    struc = parmed.Structure()
    struc.title = mol.name

    pmedatoms = []
    for atom in mol.atoms:
        pmedatm = parmed.Atom(atomic_number=atom.atomic_number,
                              name=atom.name,
                              mass=atom.mass.value_in(u.dalton),
                              number=utils.if_not_none(atom.pdbindex, -1))
        pmedatm.xx, pmedatm.xy, pmedatm.xz = atom.position.value_in(u.angstrom)
        pmedatoms.append(pmedatm)
        struc.add_atom(pmedatm,
                       resname=utils.if_not_none(atom.residue.resname, 'UNL'),
                       resnum=utils.if_not_none(atom.residue.pdbindex, -1),
                       chain=utils.if_not_none(atom.chain.name, ''))

    for bond in mol.bonds:
        struc.bonds.append(parmed.Bond(pmedatoms[bond.a1.index],
                                       pmedatoms[bond.a2.index],
                                       order=bond.order))
    return struc


def _reassign_chains(f, mol):
    """ Change chain ID assignments to the mmCIF standard (parmed uses author assignments)

    If the required fields don't exist, a copy of the molecule is returned unchanged.

    Args:
        f (file): mmcif file/stream
        mol (moldesign.Molecule): molecule with default parmed assignemnts

    Returns:
        moldesign.Molecule: new molecule with reassigned chains
    """
    data = mdt.interfaces.biopython_interface.get_mmcif_data(f)
    f.seek(0)
    try:
        newchain_names = set(data['_pdbx_poly_seq_scheme.asym_id']+
                             data['_pdbx_nonpoly_scheme.asym_id'])
    except KeyError:
        return mol.copy(name=mol.name)
    newchains = {name: mdt.Chain(name) for name in newchain_names}

    residue_iterator = itertools.chain(
            zip(data['_pdbx_poly_seq_scheme.mon_id'],
                data['_pdbx_poly_seq_scheme.pdb_seq_num'],
                data['_pdbx_poly_seq_scheme.pdb_strand_id'],
                data['_pdbx_poly_seq_scheme.asym_id']),

            zip(data['_pdbx_nonpoly_scheme.mon_id'],
                data['_pdbx_nonpoly_scheme.pdb_seq_num'],
                data['_pdbx_nonpoly_scheme.pdb_strand_id'],
                data['_pdbx_nonpoly_scheme.asym_id']))

    reschains = {(rname, ridx, rchain): newchains[chainid]
                 for rname, ridx, rchain, chainid in residue_iterator}

    for residue in mol.residues:
        newchain = reschains[residue.resname, str(residue.pdbindex), residue.chain.name]

        residue.chain = newchain

    return mdt.Molecule(mol.atoms,
                        name=mol.name, metadata=mol.metadata)

