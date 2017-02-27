""" Test various functionalities around data structure consistency
"""
import pickle

import moldesign.exceptions
import moldesign.molecules.atomcollections
import moldesign.utils.classes

from .test_objects import *


def test_h2_protected_atom_arrays(h2):
    atom1, atom2 = h2.atoms
    with pytest.raises(TypeError):
        atom1.position = 'airplane'
    with pytest.raises(u.DimensionalityError):
        atom1.momentum = 3.0 * u.femtoseconds
    with pytest.raises(ValueError):
        atom2.position = np.array([1.0, 2.0]) * u.angstrom
    atom2.position = 4.0 * u.ureg.meters
    assert atom2.x == atom2.y == atom2.z == 4.0 * u.ureg.meters


def test_h2_hierarchy(h2):
    assert len(h2.residues) == 1
    assert len(h2.chains) == 1
    chain = h2.chains.__iter__().next()
    res = h2.residues.__iter__().next()
    atom1, atom2 = h2.atoms
    assert h2 == atom1.molecule == atom2.molecule == chain.molecule == res.molecule
    assert chain == atom1.chain == atom2.chain
    assert res == atom1.residue == atom2.residue


def test_h2_array_link(h2):
    atom1, atom2 = h2.atoms
    atom2.momentum[1] = 3.0*u.default.momentum
    h2.positions[0, 1] = 0.1*u.angstrom
    assert atom1.index == 0 and atom2.index == 1
    assert atom1.y == 0.1*u.angstrom
    assert h2.momenta[1, 1] == 3.0*u.default.momentum
    assert h2.atoms[1].py == 3.0*u.default.momentum


def test_h2_harmonic_oscillator(h2_harmonic):
    mol = h2_harmonic
    atoms = h2_harmonic.atoms
    atoms[0].x = -1.0*u.angstrom
    atoms[1].x = 0.0*u.angstrom
    atoms[1].y = 0.3 * u.angstrom
    e1 = mol.calc_potential_energy()
    f1 = mol.forces[0,0]
    atoms[0].x = 1.0*u.angstrom
    f2 = mol.calc_forces()[0,0]
    e2 = mol.potential_energy

    assert e1 == e2 == 0.5*u.kcalpermol
    force_units = u.kcalpermol/u.angstrom
    assert abs(f1 + f2).value_in(force_units) < 1e-14
    assert abs(f1 - 1.0*force_units).value_in(force_units) < 1e-14


def test_h2_cache_flush(h2_harmonic):
    h2 = h2_harmonic
    pe = h2.calc_potential_energy()
    f = h2.forces
    h2.atoms[0].x += 0.1285*u.ang
    pe2 = h2.calc_potential_energy()
    f2 = h2.forces
    assert pe != pe2
    assert not np.array_equal(f, f2)


def test_h2_not_calculated_yet(h2_harmonic):
    h2_harmonic.calculate()
    h2_harmonic.atoms[1].x += 0.3*u.ang
    with pytest.raises(moldesign.exceptions.NotCalculatedError):
        h2_harmonic.forces
    with pytest.raises(moldesign.exceptions.NotCalculatedError):
        h2_harmonic.potential_energy


def h2_properties_raises_not_calculated_yet(h2_harmonic):
    h2_harmonic.calculate()
    h2_harmonic.atoms[1].x += 0.3*u.ang
    with pytest.raises(moldesign.exceptions.NotCalculatedError):
        h2_harmonic.properties.forces
    with pytest.raises(moldesign.exceptions.NotCalculatedError):
        h2_harmonic.properties.potential_energy


def test_h2_traj_energies(h2_trajectory):
    traj = h2_trajectory
    assert (np.abs(traj.positions[0,0]) <= 1.0 * u.angstrom).all()
    assert (traj.potential_energy <= 4.0 * u.kcalpermol).all()



@pytest.mark.parametrize('molkey', registered_types['molecule'])
def test_molecule_atom_hierarchy(molkey, request):
    mol = request.getfuncargvalue(molkey)
    all_residues = set(mol.residues)
    all_chains = set(mol.chains)

    residues_from_atoms = set()
    chains_from_atoms = set()
    for iatom, atom in enumerate(mol.atoms):
        assert atom.index == iatom
        assert atom.molecule is mol
        assert atom.residue in all_residues
        assert atom.chain in all_chains
        residues_from_atoms.add(atom.residue)
        chains_from_atoms.add(atom.chain)
    assert residues_from_atoms == all_residues
    assert chains_from_atoms == all_chains


@pytest.mark.parametrize('molkey', registered_types['molecule'])
def test_molecule_residue_hierarchy(molkey, request):
    mol = request.getfuncargvalue(molkey)
    all_atoms = set(mol.atoms)
    all_residues = set(mol.residues)
    all_chains = set(mol.chains)

    assert len(all_residues) == len(mol.residues)
    atoms_in_residues = set()
    chains_from_residues = set()
    assert len(all_atoms) == mol.num_atoms == len(mol.atoms)

    for residue in mol.residues:
        assert residue.molecule is mol
        chains_from_residues.add(residue.chain)
        for atom in residue.iteratoms():
            assert atom not in atoms_in_residues, 'Atom in more than 1 residue'
            atoms_in_residues.add(atom)
            assert atom in all_atoms
            assert atom.residue is residue
            assert residue.chain in all_chains

    assert chains_from_residues == all_chains
    assert atoms_in_residues == all_atoms


@pytest.mark.parametrize('molkey', registered_types['molecule'])
def test_molecule_chain_hierarchy(molkey, request):
    mol = request.getfuncargvalue(molkey)
    all_residues = set(mol.residues)
    all_chains = set(mol.chains)

    assert len(all_chains) == len(mol.chains)
    residues_from_chains = set()

    for chain in mol.chains:
        assert chain.molecule is mol
        for residue in chain:
            assert residue not in residues_from_chains, 'Residue in more than 1 chain'
            residues_from_chains.add(residue)
            assert residue in all_residues
            assert residue.molecule is mol
            assert residue.chain is chain

    assert residues_from_chains == all_residues


@pytest.mark.parametrize('molkey', registered_types['molecule'])
def test_molecule_bonds(molkey, request):
    mol = request.getfuncargvalue(molkey)
    all_atoms = set(mol.atoms)
    for atom in all_atoms:
        if atom not in mol.bond_graph:
            assert not atom.bond_graph
            continue
        bonds = mol.bond_graph[atom]
        assert atom.bond_graph == bonds
        for nbr in bonds:
            assert nbr.bond_graph[atom] == atom.bond_graph[nbr]


@pytest.mark.parametrize('molkey', registered_types['molecule'])
def test_molecule_types(molkey, request):
    mol = request.getfuncargvalue(molkey)
    assert issubclass(type(mol.atoms), moldesign.molecules.atomcollections.AtomList)
    for atom in mol.atoms:
        assert issubclass(type(atom), mdt.Atom)
    for residue in mol.residues:
        assert issubclass(type(residue), mdt.Residue)


@pytest.mark.parametrize('objkey', all_objects)
def test_pickling(objkey, request):
    obj = request.getfuncargvalue(objkey)
    for iprotocol in (0,1,2):
        x = pickle.dumps(obj, protocol=iprotocol)
        y = pickle.loads(x)
        assert type(y) == type(obj)


@pytest.mark.parametrize('objkey', registered_types['equality'])
def test_pickled_equality(objkey, request):
    obj = request.getfuncargvalue(objkey)

    for iprotocol in (0,1,2):
        x = pickle.dumps(obj, protocol=iprotocol)
        y = pickle.loads(x)
        assert type(y) == type(obj)
        try:
            assert y == obj
        except ValueError:
            assert (y == obj).all()


