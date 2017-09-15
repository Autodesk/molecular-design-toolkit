from builtins import zip

import itertools
import numpy as np
import pytest

import moldesign as mdt
from moldesign import units as u

from .molecule_fixtures import *
from .object_fixtures import *


__PYTEST_MARK__ = 'internal'  # mark all tests in this module with this label (see ./conftest.py)


def test_carbon_copy(carbon_copy, carbon_atom):
    atom = carbon_copy
    assert atom.symbol == carbon_atom.symbol
    assert atom.mass == carbon_atom.mass
    assert atom.bond_graph == {}


@pytest.mark.screening
def test_copy_breaks_link(h2):
    h2copy = mdt.Molecule(h2)
    h2.atoms[0].y = 4.0 * u.bohr
    assert h2copy.atoms[0].y == 0.0 * u.angstrom
    np.testing.assert_almost_equal(h2.positions[0,1].value_in(u.bohr),
                                   4.0, 7)
    assert h2copy.positions[0, 1] == 0.0 * u.bohr

    h2copy.momenta[1, 0] = 2.0 * u.default.momentum
    np.testing.assert_almost_equal(h2copy.atoms[1].px.value_in(u.default.momentum),
                                   2.0, 7)
    assert h2.momenta[1, 0] == 0.0 * u.default.momentum
    assert h2.atoms[1].px == 0.0 * u.default.momentum


def test_h2_harmonic_copy(h2_harmonic_copy, h2_harmonic):
    mol = h2_harmonic_copy
    assert mol.energy_model is mol.integrator is None
    for name in 'num_atoms num_residues positions momenta'.split():
        try:
            assert getattr(h2_harmonic_copy, name) == getattr(h2_harmonic, name)
        except ValueError:
            assert (getattr(h2_harmonic_copy, name) == getattr(h2_harmonic, name)).all()

    assert mol.atoms[0].bond_graph[mol.atoms[1]] == 1
    assert mol.atoms[1].bond_graph[mol.atoms[0]] == 1


@pytest.mark.screening
def test_copying_doesnt_corrupt_original_h2_harmonic(h2_harmonic):
    mol = h2_harmonic
    integ = mol.integrator
    model = mol.energy_model
    residue = mol.residues[0]
    chain = list(mol.chains)[0]
    m2 = mdt.Molecule(mol)
    assert integ is mol.integrator
    assert model is mol.energy_model
    assert len(mol.chains) == 1
    assert len(mol.residues) == 1
    assert residue == mol.residues[0]
    assert chain == list(mol.chains)[0]


def test_atoms_copied_from_h2_harmonic(copy_atoms_from_h2_harmonic, h2_harmonic):
    atoms = copy_atoms_from_h2_harmonic
    atom = atoms[0]
    assert atom.molecule is None
    assert atom.residue is not h2_harmonic.atoms[0].residue
    assert atom.chain is not h2_harmonic.atoms[0].chain
    assert atoms[0].residue is atoms[1].residue
    assert atoms[0].chain is atoms[1].chain


@pytest.mark.parametrize('fixture_key',
                         registered_types['molecule'] + registered_types['submolecule'])
def test_copy_atoms(fixture_key, request):
    oldobj = request.getfixturevalue(fixture_key)
    oldatoms = oldobj.atoms
    atoms = oldatoms.copy_atoms()
    assert isinstance(atoms, mdt.AtomList)
    assert_atom_copy_integrity(atoms, oldatoms)


@pytest.mark.parametrize('fixture_key',
                         registered_types['molecule'] + registered_types['submolecule'])
def test_copy_atoms_in_containers(fixture_key, request):
    oldobj = request.getfixturevalue(fixture_key)
    if isinstance(oldobj, mdt.molecules.ChildList):
        pytest.xfail("We haven't defined the behavior for copying ChildLists yet")

    newobj = oldobj.copy()
    assert newobj.__class__ is oldobj.__class__
    assert_atom_copy_integrity(newobj.atoms, oldobj.atoms)


def assert_atom_copy_integrity(atoms, oldatoms):
    newmol = atoms[0].molecule
    assert len(atoms) == len(oldatoms)
    for oldatom, atom in itertools.zip_longest(oldatoms, atoms):
        assert oldatom.name == atom.name
        assert oldatom.atnum == atom.atnum
        assert oldatom.mass == atom.mass
        assert (oldatom.position == atom.position).all()
        assert (oldatom.velocity == atom.velocity).all()

        atom.x += 1.0 * u.angstrom
        assert all((oldatom.position == atom.position) == [False, True, True])

        atom.px += 1.0 * u.angstrom * u.amu / u.fs
        assert all((oldatom.momentum == atom.momentum) == [False, True, True])

        assert oldatom.residue is not atom.residue
        assert oldatom.residue.name == atom.residue.name
        assert oldatom.chain is not atom.chain
        assert oldatom.chain.name == atom.chain.name
        assert (oldatom.molecule is atom.molecule is None) or (oldatom.molecule is not atom.molecule)

        assert atom.molecule is newmol  # might be None

        if oldatom.residue is not None:
            assert atom.residue.molecule is newmol
            assert atom in atom.residue
            assert atom.residue[atom.name] is atom
        if oldatom.chain is not None:
            assert atom.chain.molecule is newmol
            assert atom in atom.chain.atoms
            assert atom.residue in atom.chain
            assert atom.chain[atom.residue.name] is atom.residue
        if newmol is not None:
            assert atom.residue is newmol.residues[atom.residue.index]
            assert atom.chain is newmol.chains[atom.chain.index]
            assert atom is newmol.atoms[atom.index]


@pytest.mark.parametrize('fixture_key', moldesign_objects)
def test_copied_types(fixture_key, request):
    obj = request.getfixturevalue(fixture_key)
    if isinstance(obj, mdt.molecules.ChildList):
        pytest.xfail("We haven't defined the behavior for copying ChildLists yet")
    objcopy = obj.copy()
    assert obj.__class__ is objcopy.__class__

    if not (isinstance(obj, mdt.Trajectory) or isinstance(obj, mdt.Atom)):
        atomcopy = obj.copy_atoms()
        assert isinstance(atomcopy, mdt.AtomList)
        assert len(atomcopy) == obj.num_atoms


@pytest.mark.parametrize('fixture_key', registered_types['molecule'])
def test_molecule_copy(fixture_key, request):
    mol = request.getfixturevalue(fixture_key)
    newmol = mol.copy()
    assert mol.name in newmol.name
    assert mol.energy_model == newmol.energy_model
    assert mol.integrator == newmol.integrator


@pytest.mark.parametrize('fixturename', 'pdb3aid h2'.split())
def test_molecular_combination_chains(request, fixturename):
    mol = request.getfixturevalue(fixturename)
    m2 = mol.copy()
    newmol = mol.combine(m2)
    assert newmol.num_chains == 2*mol.num_chains
    assert len(set(chain for chain in newmol.chains)) == newmol.num_chains
    _assert_unique_chainids(newmol)


def _assert_unique_chainids(mol):
    chain_names = set()
    for chain in mol.chains:
        assert chain.name == chain.pdbindex
        assert chain.name not in chain_names
        chain_names.add(chain.name)


def test_chain_rename(pdb3aid):
    res1 = mdt.Molecule(pdb3aid.residues[3])
    res2 = mdt.Molecule(pdb3aid.residues[4])
    newmol = mdt.Molecule([res1, res2])
    assert newmol.num_chains == 2
    assert newmol.num_residues == 2
    assert newmol.residues[0].name == res1.residues[0].name
    assert newmol.residues[1].name == res2.residues[0].name
    assert newmol.chains[0].name == 'A'
    assert newmol.chains[1].name == 'B'


@pytest.mark.parametrize('molkey', ["cached_mol_parameterized_with_zeros",
                                    "cached_protein_with_default_amber_ff"])
def test_forcefield_copied_with_molecule(molkey, request):
    mol = request.getfixturevalue(molkey)
    m2 = mol.copy()

    assert isinstance(m2.ff, mol.ff.__class__)
    assert isinstance(m2.ff.parmed_obj, mol.ff.parmed_obj.__class__)
    assert m2.ff.parmed_obj is not mol.ff.parmed_obj

    p1 = m2.ff.parmed_obj
    p2 = m2.ff.parmed_obj
    assert m2.ff.parmed_obj.LJ_depth == mol.ff.parmed_obj.LJ_depth
    assert p1.bond_types == p2.bond_types
    assert p1.angle_types == p2.angle_types
    assert p1.dihedral_types == p2.dihedral_types
    assert p1.improper_types == p2.improper_types


def test_constraints_copied_with_molecule(mol_with_zerocharge_params):
    mol = mol_with_zerocharge_params

    mol.constrain_distance(*mol.atoms[:2])
    mol.constrain_angle(*mol.atoms[:3])
    mol.constrain_dihedral(*mol.atoms[:4])
    mol.constrain_atom(mol.atoms[0])
    mol.constrain_hbonds()

    mcpy = mol.copy()
    assert isinstance(mcpy.constraints, mol.constraints.__class__)
    assert mcpy.constraints[0].desc == 'distance'
    assert mcpy.constraints[1].desc == 'angle'
    assert mcpy.constraints[2].desc == 'dihedral'
    assert mcpy.constraints[3].desc == 'position'
    assert mcpy.constraints[4].desc == 'hbonds'

    for constraint in mcpy.constraints:
        assert constraint.mol is mcpy
        if constraint.desc != 'hbonds':
            for atom in constraint.atoms:
                assert atom.molecule is mcpy


def test_properties_copied_with_molecule(cached_h2_rhf_sto3g):
    original = cached_h2_rhf_sto3g
    assert original.potential_energy is not None  # sanity check

    mol = cached_h2_rhf_sto3g.copy()

    assert mol is not original
    assert mol.properties is not original.properties

    for prop, val in original.properties.items():
        assert prop in mol.properties
        if isinstance(val, (u.MdtQuantity, np.ndarray)) and getattr(val, 'shape', False):
                assert (mol.properties[prop] == val).all()
                assert mol.properties[prop] is not val

        elif isinstance(val, (str, int, float, mdt.molecules.AtomicProperties)):
            assert mol.properties[prop] == val

        else:  # otherwise, just make sure it's not the original
            assert mol.properties[prop] is not val


def test_wfn_copied_with_molecule(cached_h2_rhf_sto3g):
    original = cached_h2_rhf_sto3g
    assert original.wfn is not None  # sanity check

    mol = original.copy()

    assert mol.wfn is not None

    # should be completely equal
    assert (mol.wfn.aobasis.fock == original.wfn.aobasis.fock).all()
    # but different objects
    assert mol.wfn.aobasis.fock is not original.wfn.aobasis.fock


def test_wfn_copy(cached_h2_rhf_sto3g):
    original = cached_h2_rhf_sto3g
    wfn = original.wfn.copy()

    assert wfn.mol is original
    assert wfn is not original.wfn

    # should be completely equal
    assert (wfn.aobasis.fock == original.wfn.aobasis.fock).all()
    # but different objects
    assert wfn.aobasis.fock is not original.wfn.aobasis.fock