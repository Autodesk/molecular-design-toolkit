import numpy as np
import pytest

import moldesign as mdt
from moldesign import units as u

from .test_objects import *


def test_carbon_copy(carbon_copy, carbon_atom):
    atom = carbon_copy
    assert atom.symbol == carbon_atom.symbol
    assert atom.mass == carbon_atom.mass
    assert atom.bond_graph == {}


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


def test_h2_harmonic_copy_loses_simulation(h2_harmonic_copy, h2_harmonic):
    mol = h2_harmonic_copy
    assert mol.energy_model is mol.integrator is None
    for name in 'num_atoms num_residues positions momenta'.split():
        try:
            assert getattr(h2_harmonic_copy, name) == getattr(h2_harmonic, name)
        except ValueError:
            assert (getattr(h2_harmonic_copy, name) == getattr(h2_harmonic, name)).all()

    assert mol.atoms[0].bond_graph[mol.atoms[1]] == 1
    assert mol.atoms[1].bond_graph[mol.atoms[0]] == 1


def test_h2_calculation_caching(h2_harmonic):
    h2 = h2_harmonic
    h2.properties = mdt.MolecularProperties(h2)
    true_energy = h2.calc_potential_energy()
    assert 'potential_energy' in h2.properties
    assert 'forces' in h2.properties
    h2.potential_energy
    h2.forces
    h2.properties['potential_energy'] = 'banana'
    assert h2.potential_energy == h2.calc_potential_energy() == 'banana'
    props = h2.calculate()
    assert props.potential_energy == h2.potential_energy == h2.calc_potential_energy() == 'banana'
    props2 = h2.calculate(use_cache=False)
    assert props2.potential_energy == h2.potential_energy == true_energy
    assert h2.calc_potential_energy() == true_energy


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


@pytest.mark.parametrize('fixture_key', ['h2_harmonic_atoms',
                                         'ligand_3aid_atoms',
                                         'random_atoms_from_3aid'])
def test_copy_atoms(fixture_key, request):
    oldatoms = request.getfuncargvalue(fixture_key)
    atoms = oldatoms.copy()
    assert isinstance(atoms, mdt.AtomList)
    for old, new in zip(oldatoms, atoms):
        assert old.residue is not new.residue
        assert old.residue.name == new.residue.name
        assert old.chain is not new.chain
        assert old.chain.name == new.chain.name
        assert new.molecule is None


@pytest.mark.parametrize('fixture_key', moldesign_objects)
def test_copied_types(fixture_key, request):
    obj = request.getfuncargvalue(fixture_key)
    objcopy = obj.copy()
    assert obj.__class__ is objcopy.__class__

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


