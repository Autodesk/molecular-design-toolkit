""" Tests Molecule instances and other base classes they depend on
"""

import pickle
import random

import numpy as np
import pytest

import moldesign as mdt
import moldesign.exceptions
import moldesign.molecules.atomcollections
import moldesign.utils.classes
from moldesign import units as u

registered_types = {}


def typedfixture(*types, **kwargs):
    """This is a decorator that lets us associate fixtures with one or more arbitrary types.
    We'll later use this type to determine what tests to run on the result"""

    def fixture_wrapper(func):
        for t in types:
            registered_types.setdefault(t, []).append(func.__name__)
        return pytest.fixture(**kwargs)(func)

    return fixture_wrapper


######################################
# Test the basic data structures
@typedfixture('object')
def dotdict():
    dd = moldesign.utils.classes.DotDict(a='b',
                                         c=3,
                                         d='e')
    return dd


@typedfixture('object')
def dictlike():
    dd = moldesign.utils.classes.DictLike(a='b',
                                          c=3,
                                          d='e')
    return dd


def test_dotdict(dotdict):
    dd = dotdict
    assert dd.a == 'b'
    assert dd.d == 'e'


@pytest.mark.parametrize('objkey', 'dictlike dotdict'.split())
def test_dictlike(objkey, request):
    dd = request.getfuncargvalue(objkey)
    assert set(dd.keys()) == {'a', 'c', 'd'}
    assert set(dd.values()) == {'b', 3, 'e'}


# Some objects with units
@typedfixture('object')
def list_of_units(): return [1.0 * u.angstrom, 1.0 * u.nm, 1.0 * u.a0]


@typedfixture('object', 'equality')
def simple_unit_array(): return np.array([1.0, -2.0, 3.5]) * u.angstrom


@typedfixture('object', 'equality')
def unit_number(): return 391.23948 * u.ureg.kg * u.ang / u.alpha


######################################
# Test underlying elements
@typedfixture('submolecule')
def carbon_atom():
    atom1 = mdt.Atom('C')
    return atom1


def test_carbon_atom(carbon_atom):
    assert carbon_atom.symbol == 'C'
    assert carbon_atom.mass == 12.0 * u.amu


@typedfixture('submolecule')
def carbon_copy(carbon_atom):
    atoms = carbon_atom.copy()
    return atoms


def test_carbon_copy(carbon_copy, carbon_atom):
    atom = carbon_copy
    assert atom.symbol == carbon_atom.symbol
    assert atom.mass == carbon_atom.mass
    assert atom.bond_graph == {}


######################################
# Tests around hydrogen
@typedfixture('molecule')
def h2():
    atom1 = mdt.Atom('H')
    atom1.x = 0.5 * u.angstrom
    atom2 = mdt.Atom(atnum=1)
    atom2.position = [-0.5, 0.0, 0.0] * u.angstrom
    h2 = mdt.Molecule([atom1, atom2], name='h2')
    atom1.bond_to(atom2, 1)
    return h2


@typedfixture('molecule')
def h2_harmonic():
    mol = h2()
    SPRING_CONSTANT = 1.0 * u.kcalpermol / (u.angstrom ** 2)
    model = moldesign.models.HarmonicOscillator(k=SPRING_CONSTANT)
    integrator = moldesign.integrators.VelocityVerlet(timestep=0.5*u.fs, frame_interval=30)
    mol.set_energy_model(model)
    mol.set_integrator(integrator)
    return mol


@typedfixture('submolecule')
def h2_trajectory(h2_harmonic):
    mol = h2_harmonic
    mol.atoms[0].x = 1.0 * u.angstrom
    mol.momenta *= 0.0
    traj = mol.run(500)
    return traj


def test_h2_trajectory(h2_trajectory):
    """
    Check if the trajectory is the sine wave that we expect
    """
    traj = h2_trajectory
    mol = traj.mol
    k = mol.energy_model.params.k
    period = 2*u.pi*np.sqrt(mol.atoms[0].mass/k)
    for frame in traj.frames:
        period_progress = (frame.time % period) / period
        if period_progress < 0.1 or period_progress > 0.9:
            # check for expected peaks of sine wave
            assert frame.positions[0, 0] > 0.1 * u.angstrom
        elif 0.4 < period_progress < 0.6:
            # check for expected troughs of sine wave
            assert frame.positions[0, 0] < -0.1 * u.angstrom


@typedfixture('molecule')
def h2_traj_tempmol(h2_trajectory):
    return h2_trajectory._tempmol


@typedfixture('molecule')
def h2_harmonic_copy(h2_harmonic):
    return mdt.Molecule(h2_harmonic)


def test_h2_positions(h2):
    atom1, atom2 = h2.atoms
    assert (atom1.position == np.array([0.5, 0.0, 0.0]) * u.angstrom).all()
    assert atom2.x == -0.5 * u.angstrom
    assert atom1.distance(atom2) == 1.0 * u.angstrom


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
    h2.properties = moldesign.molecules.molecule.MolecularProperties(h2)
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


@typedfixture('submolecule')
def copy_atoms_from_h2_harmonic(h2_harmonic):
    atoms = h2_harmonic.atoms.copy()
    return atoms

@typedfixture('molecule')
def h2_harmonic_thats_been_copied(h2_harmonic):
    temp = mdt.Molecule(h2_harmonic)
    return h2_harmonic


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


@typedfixture('submolecule')
def h2_harmonic_atoms(h2_harmonic):
    return h2_harmonic.atoms


def test_h2_traj_energies(h2_trajectory):
    traj = h2_trajectory
    assert (np.abs(traj.positions[0,0]) <= 1.0 * u.angstrom).all()
    assert (traj.potential_energy <= 4.0 * u.kcalpermol).all()


def test_h2(h2):
    mol = h2
    assert mol.num_atoms == 2
    assert mol.atoms[0].symbol == mol.atoms[1].symbol == 'H'


######################################
# Tests around PDB ID 3AID
@typedfixture('molecule', scope='session')
def pdb3aid():
    #FPATH = '../notebooks/data/3AID.pdb'
    #assert os.path.exists(FPATH)
    mol = mdt.from_pdb('3AID')
    return mol


@typedfixture('submolecule')
def ligand_residue_3aid(pdb3aid):
    unknown = pdb3aid.chains['A'](type='unknown')
    assert len(unknown) == 1
    return unknown[0]


@typedfixture('submolecule')
def ligand_3aid_atoms(ligand_residue_3aid):
    return ligand_residue_3aid.atoms


@typedfixture('molecule')
def ligand3aid(ligand_residue_3aid):
    newmol = mdt.Molecule(ligand_residue_3aid)
    return newmol


@pytest.fixture
def random_atoms_from_3aid(pdb3aid):
    atoms = moldesign.molecules.atomcollections.AtomList(random.sample(pdb3aid.atoms, 10))
    return atoms


def test_3aid(pdb3aid):
    mol = pdb3aid
    assert len(mol.chains) == 2


def test_3aid_ligand_search(pdb3aid):
    mol = pdb3aid
    unknown = mol.chains['A'](type='unknown')
    proteina = mol.chains['A'](type='protein')
    proteinb = mol.chains['B'](type='protein')
    assert len(unknown) == 1
    assert len(proteina) == len(proteinb) == 99


def test_ligand3aid(ligand3aid):
    mol = ligand3aid
    assert len(mol.chains) == 1
    assert len(mol.residues) == 1


def test_nucleic_build(nucleic):
    mol = nucleic
    assert mol.num_chains == 2
    assert mol.num_residues == 8
    assert mol.chains[0] is mol.chains['A']
    assert mol.chains[1] is mol.chains['B']
    assert len(mol.chains[0].residues) == len(mol.chains[1].residues) == 4



######################################
# Tests around a piece of DNA
@typedfixture('molecule', scope='session')
def nucleic():
    # ACTG.pdb contains a molecule generated using mdt.build_dna('ACTG')
    mol = mdt.read('data/ACTG.pdb')
    return mol


######################################
# Consistency tests for groups of objects
moldesign_objects = registered_types['molecule'] + registered_types['submolecule']
all_objects = registered_types['object'] + moldesign_objects


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
    x = pickle.dumps(obj, protocol=2)
    y = pickle.loads(x)
    assert type(y) == type(obj)


@pytest.mark.parametrize('objkey', registered_types['equality'])
def test_pickled_equality(objkey, request):
    obj = request.getfuncargvalue(objkey)
    x = pickle.dumps(obj, protocol=2)
    y = pickle.loads(x)
    assert type(y) == type(obj)
    try:
        assert y == obj
    except ValueError:
        assert (y == obj).all()


@pytest.mark.parametrize('fixture_key', ['h2_harmonic_atoms',
                                         'ligand_3aid_atoms',
                                         'random_atoms_from_3aid'])
def test_copy_atoms(fixture_key, request):
    oldatoms = request.getfuncargvalue(fixture_key)
    atoms = oldatoms.copy()
    for old, new in zip(oldatoms, atoms):
        assert old.residue is not new.residue
        assert old.residue.name == new.residue.name
        assert old.chain is not new.chain
        assert old.chain.name == new.chain.name

