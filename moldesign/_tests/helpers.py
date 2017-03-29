import os

import numpy as np

import moldesign as mdt
from moldesign import units as u

DEFSTEP = 0.000005*u.angstrom


def get_data_path(f):
    """
    Returns path to a file in the ``moldesign/_tests/data/`` directory.
    """
    return os.path.join(mdt.data.PACKAGEPATH, '_tests', 'data', f)


def num_grad(mol, fn, step=DEFSTEP, atoms=None, fnargs=None, fnkwargs=None):
    """ Calculate the finite-difference gradient of a function
    """
    grad = None
    origpos = mol.positions.copy()
    if fnargs is None:
        fnargs = tuple()
    if fnkwargs is None:
        fnkwargs = dict()

    if atoms is None:
        atoms = mol.atoms

    for iatom, atom in enumerate(atoms):
        for idim in xrange(3):
            atom.position[idim] += step
            vplus = fn(*fnargs, **fnkwargs)
            atom.position[idim] -= 2.0 * step
            vminus = fn(*fnargs, **fnkwargs)
            mol.positions = origpos  # reset positions

            if grad is None:
                grad = np.zeros((len(atoms), 3)) * vplus.units/mol.positions.units
            grad[iatom, idim] = (vplus - vminus) / (2.0*step)

    return grad


def _make_mol_with_n_hydrogens(n):
    return mdt.Molecule([mdt.Atom('H') for i in xrange(n)])


class ZeroEnergy(mdt.models.base.EnergyModelBase):
    """ All 0, all the time
    """
    def prep(self):
        pass

    def calculate(self):
        return dict(potential_energy=0.0*u.default.energy,
                    forces=np.zeros(self.mol.positions.shape)*u.default.force)


def check_topology_consistency(mol):
    """ General helper to make sure the data hasn't gotten corrupted. Use this liberally during
    testing!
    """
    all_residues = set(mol.residues)
    populated_residues = set([residue for residue in mol.residues if residue.num_atoms>0])
    all_chains = set(mol.chains)

    populated_chains = set([chain for chain in all_chains if chain.num_atoms>0])
    found_chains = set(atom.chain for atom in mol.atoms)
    found_residues = set(atom.residue for atom in mol.atoms)

    assert len(all_chains) == mol.num_chains
    assert found_chains == populated_chains
    assert len(all_residues) == mol.num_residues
    assert found_residues == populated_residues

    residues_from_chains = set()
    for ichain, chain in enumerate(mol.chains):
        assert chain.index == ichain
        assert mol.chains[chain.index] is chain
        if chain.name is not None:
            assert mol.chains[chain.name] is chain
        assert chain.molecule is mol

        for residue in chain:
            assert residue not in residues_from_chains, 'Residue in more than 1 chain'
            residues_from_chains.add(residue)
            assert residue in all_residues
            assert residue.molecule is mol
            assert residue.chain is chain
            assert mol.residues[residue.index] is residue
            if residue.name is not None:
                assert chain[residue.name] is residue

            prevatom = None
            for atom in residue:
                assert atom.residue is residue
                assert atom.chain is chain
                assert atom.molecule is mol
                assert mol.atoms[atom.index] is atom
                assert residue[atom.name] is atom
                if prevatom is not None:
                    assert (prevatom.index, prevatom.pdbindex) < (atom.index, atom.pdbindex)
                prevatom = atom

    assert residues_from_chains == all_residues

    check_bond_consistency(mol)


def check_bond_consistency(mol):
    """ Checks the internal structure of the bond graph is consistent
    """
    found_bonds = 0
    for bond in mol.bonds:
        found_bonds += 1
        assert mol.bonds[bond.a2][bond.a1] == bond
        assert mol.bonds[bond.a1][bond.a2] == bond
        assert bond.a1 in mol
        assert bond.a2 in mol
    assert found_bonds == mol.num_bonds

    num_atombonds = 0
    for atom in mol.atoms:

        bondset = set()
        for bond in atom.bonds:
            bondset.add(bond)
            nbr = bond.partner(atom)
            assert nbr in mol
            assert nbr in atom.bonds
            num_atombonds += 1

        assert len(bondset) == atom.num_bonds
        num_nbr = 0
        for nbr in atom.bonds.atoms:
            num_nbr += 1
            b = mdt.Bond(atom, nbr)
            assert b in bondset
            assert b == atom.bonds[nbr]
            assert b == mol.bonds[atom][nbr]
            assert b == mol.bonds[nbr][atom]

        assert num_nbr == atom.num_bonds
        assert len(atom.bonds) == len(mol.bonds[atom]) == atom.num_bonds

    assert num_atombonds == 2 * mol.num_bonds
