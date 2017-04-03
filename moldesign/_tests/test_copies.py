import numpy as np
import pytest

import moldesign as mdt
from moldesign import units as u

from .object_fixtures import *


def test_carbon_copy(carbon_copy, carbon_atom):
    atom = carbon_copy
    assert atom.symbol == carbon_atom.symbol
    assert atom.mass == carbon_atom.mass
    assert atom.num_bonds == 0


def test_copy_breaks_link(h2):
    h2copy = mdt.Molecule(h2)
    h2.atoms[0].y = 4.0 * u.bohr
    assert h2copy.atoms[0].y == 0.0 * u.angstrom
    np.testing.assert_almost_equal(h2.positions[0,1].value_in(u.bohr),
                                   4.0, 7)
    assert h2copy.positions[0, 1] == 0.0 * u.bohr

    assert h2copy.atoms[0] in h2copy.atoms[1].bonds

    h2copy.momenta[1, 0] = 2.0 * u.default.momentum
    np.testing.assert_almost_equal(h2copy.atoms[1].px.value_in(u.default.momentum),
                                   2.0, 7)
    assert h2.momenta[1, 0] == 0.0 * u.default.momentum
    assert h2.atoms[1].px == 0.0 * u.default.momentum


def test_h2_harmonic_copy_loses_simulation(h2_harmonic_copy, h2_harmonic):
    mol = h2_harmonic_copy
    assert mol.energy_model is mol.integrator is None
    for name in 'num_atoms num_residues positions momenta'.split():
        try:
            assert getattr(h2_harmonic_copy, name) == getattr(h2_harmonic, name)
        except ValueError:
            assert (getattr(h2_harmonic_copy, name) == getattr(h2_harmonic, name)).all()

    assert mol.atoms[0].bonds[mol.atoms[1]].order == 1
    assert mol.atoms[1].bonds[mol.atoms[0]].order == 1


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
    assert atoms[0] in atoms[1].bonds
    assert atoms[1] in atoms[0].bonds


@pytest.mark.parametrize('fixture_key',
                         registered_types['molecule'] + registered_types['submolecule'])
def test_copy_atoms(fixture_key, request):
    oldobj = request.getfuncargvalue(fixture_key)
    oldatoms = oldobj.atoms
    atoms = oldatoms.copy_atoms()
    assert isinstance(atoms, mdt.AtomList)
    _test_copy_integrity(atoms, oldatoms)


@pytest.mark.parametrize('fixture_key',
                         registered_types['molecule'] + registered_types['submolecule'])
def test_copy_atoms_in_containers(fixture_key, request):
    oldobj = request.getfuncargvalue(fixture_key)
    if isinstance(oldobj, mdt.molecules.ChildList):
        pytest.xfail("We haven't defined the behavior for copying ChildLists yet")

    newobj = oldobj.copy()
    assert newobj.__class__ is oldobj.__class__
    _test_copy_integrity(newobj.atoms, oldobj.atoms)


def _test_copy_integrity(atoms, oldatoms):
    for old, new in zip(oldatoms, atoms):
        assert (old.position == new.position).all()
        assert (old.velocity == new.velocity).all()

        new.x += 1.0 * u.angstrom
        assert all((old.position == new.position) == [False, True, True])

        new.px += 1.0 * u.angstrom * u.amu / u.fs
        assert all((old.momentum == new.momentum) == [False, True, True])

        assert old.residue is not new.residue
        assert old.residue.name == new.residue.name
        assert old.chain is not new.chain
        assert old.chain.name == new.chain.name
        assert (old.molecule is new.molecule is None) or (old.molecule is not new.molecule)


@pytest.mark.parametrize('fixture_key', moldesign_objects)
def test_copied_types(fixture_key, request):
    obj = request.getfuncargvalue(fixture_key)
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
    mol = request.getfuncargvalue(fixture_key)
    newmol = mol.copy()
    assert mol.name in newmol.name
    assert mol.energy_model == newmol.energy_model
    assert mol.integrator == newmol.integrator


def test_molecular_combination(pdb3aid):
    m2 = pdb3aid.copy()
    newmol = pdb3aid.combine(m2)
    assert newmol.num_chains == 2 * pdb3aid.num_chains
    assert newmol.num_residues == 2 * pdb3aid.num_residues
    assert newmol.num_atoms == 2 * pdb3aid.num_atoms
    assert newmol.num_bonds == 2 * pdb3aid.num_bonds
    assert len(set(chain for chain in newmol.chains)) == newmol.num_chains


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

