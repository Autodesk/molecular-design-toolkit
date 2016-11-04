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

import itertools
import parmed

import moldesign as mdt
from moldesign import units as u


def exports(o):
    __all__.append(o.__name__)
    return o
__all__ = []


def parse_mmcif(f, reassign_chains=True):
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


def parse_pdb(f):
    """Parse an mmCIF file (using the parmEd parser) and return a molecule

    Args:
        f (file): file-like object containing the mmCIF file

    Returns:
        moldesign.Molecule: parsed molecule
    """
    parmedmol = parmed.read_PDB(f)
    mol = parmed_to_mdt(parmedmol)
    return mol


@exports
def parmed_to_mdt(pmdmol):
    """ Convert parmed object to MDT object

    Args:
        pmdmol (parmed.Structure): parmed structure to convert

    Returns:
        mdt.Molecule: converted molecule
    """
    atoms = collections.OrderedDict()
    residues = {}
    chains = {}
    for patm in pmdmol.atoms:
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
                        mass=patm.mass * u.dalton)
        atom.position = [patm.xx, patm.xy, patm.xz]*u.angstrom

        atom.residue = residue
        residue.add(atom)
        atom.chain = chain
        atoms[patm] = atom

    for pbnd in pmdmol.bonds:
        atoms[pbnd.atom1].bond_to(atoms[pbnd.atom2], int(pbnd.order))

    mol = mdt.Molecule(atoms.values())
    mol.description = pmdmol.title
    return mol


def _parmed_to_ff(topo, atom_map):
    """ Create an MDT FFParameters object from a ParmEd topology

    Args:
        topo (parmed.Structure): ParmEd structure (with FF terms)
        atom_map (Mapping[parmed.Atom, moldesign.Atom]): mapping between MDT and ParmEd atoms

    Returns:
        moldesign.forcefields.FFParameters: parameters in MDT format
    """
    bonds = [mdt.forcefields.HarmonicBondTerm(atom_map[bond.a1],
                                              atom_map[bond.a2],
                                              bond.type.k*u.kcalpermol/u.angstrom ** 2,
                                              bond.type.req*u.angstrom)
             for bond in topo.bonds]

    angles = [mdt.forcefields.HarmonicAngleTerm(atom_map[angle.a1],
                                                atom_map[angle.a2],
                                                atom_map[angle.a3],
                                                angle.type.k*u.kcalpermol/u.radian ** 2,
                                                angle.type.theta_eq*u.degrees)
              for angle in topo.angles]

    dihedrals = [mdt.forcefields.PeriodicTorsionTerm(atom_map[dihedral.a1],
                                                     atom_map[dihedral.a2],
                                                     atom_map[dihedral.a3],
                                                     atom_map[dihedral.a4],
                                                     dihedral.type.per,
                                                     dihedral.type.phi_k*u.kcalpermol,
                                                     dihedral.type.phase*u.degrees)
                 for dihedral in topo.dihedrals]


def _reassign_chains(f, mol):
    """ Change chain ID assignments to the mmCIF standard (parmed uses author assignments)

    Args:
        f (file): mmcif file/stream
        mol (moldesign.Molecule): molecule with default parmed assignemnts

    Returns:
        moldesign.Molecule: new molecule with reassigned chains
    """
    data = mdt.interfaces.biopython_interface.get_mmcif_data(f)
    f.seek(0)
    newchain_names = set(data['_pdbx_poly_seq_scheme.asym_id']+
                         data['_pdbx_nonpoly_scheme.asym_id'])
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

        for atom in residue.atoms:
            atom.chain = newchain
        residue.chain = newchain

    return mdt.Molecule(mol.atoms,
                        name=mol.name, description=mol.description)

