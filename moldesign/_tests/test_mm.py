import random

import pytest
import numpy as np

import moldesign as mdt
from moldesign import units as u

from . import helpers

registered_types = {}


def typedfixture(*types, **kwargs):
    """This is a decorator that lets us associate fixtures with one or more arbitrary types.
    We'll later use this type to determine what tests to run on the result"""

    def fixture_wrapper(func):
        for t in types:
            registered_types.setdefault(t, []).append(func.__name__)
        return pytest.fixture(**kwargs)(func)

    return fixture_wrapper


@pytest.fixture
def small_molecule():
    mol = mdt.from_smiles('CNCOS(=O)C')
    mol.positions += 0.001*np.random.random(mol.positions.shape)*u.angstrom  # move out of minimum
    return mol


@typedfixture('hasmodel')
def parameterize_zeros(small_molecule):
    params = mdt.parameterize(small_molecule, charges='zero')
    mol = mdt.assign_forcefield(small_molecule, parameters=params)
    mol.set_energy_model(mdt.models.ForceField)
    return mol


@typedfixture('hasmodel')
def parameterize_am1bcc(small_molecule):
    params = mdt.parameterize(small_molecule, charges='am1-bcc', ffname='gaff')
    mol = mdt.assign_forcefield(small_molecule, parameters=params)
    mol.set_energy_model(mdt.models.ForceField)
    return mol


@typedfixture('hasmodel')
def openbabel_mmff94(small_molecule):
    small_molecule.set_energy_model(mdt.models.OpenBabelPotential, forcefield='mmff94')
    return small_molecule


@typedfixture('hasmodel')
def openbabel_mmff94s(small_molecule):
    small_molecule.set_energy_model(mdt.models.OpenBabelPotential, forcefield='mmff94s')
    return small_molecule


# This test (along with the uff energy model) is disabled because it does not appear to return a
# gradient that's consistent with the energy surface
#@typedfixture('hasmodel')
#def openbabel_uff(small_molecule):
#    small_molecule.set_energy_model(mdt.models.OpenBabelPotential, forcefield='uff')
#    return small_molecule


@typedfixture('hasmodel')
def openbabel_ghemical(small_molecule):
    small_molecule.set_energy_model(mdt.models.OpenBabelPotential, forcefield='ghemical')
    return small_molecule


@typedfixture('hasmodel')
def protein_default_amber_forcefield():
    mol = mdt.from_pdb('1YU8')
    newmol = mdt.assign_forcefield(mol)
    newmol.set_energy_model(mdt.models.ForceField)
    return newmol


@typedfixture('hasmodel')
def gaff_model_gasteiger(small_molecule):
    small_molecule.set_energy_model(mdt.models.GAFF, partial_charges='gasteiger')
    return small_molecule


@pytest.mark.parametrize('objkey', registered_types['hasmodel'])
def test_forces_and_energy_were_calculated(objkey, request):
    mol = request.getfuncargvalue(objkey)
    energy = mol.calculate_potential_energy()
    forces = mol.calculate_forces()
    assert forces.shape == mol.positions.shape


@pytest.mark.skipif(mdt.interfaces.openmm.force_remote,
                    reason="Numerical derivatives need to be parallelized, "
                           "otherwise this takes too long")
@pytest.mark.parametrize('objkey', registered_types['hasmodel'])
def test_analytical_vs_numerical_forces(objkey, request):
    mol = request.getfuncargvalue(objkey)

    if mol.num_atoms > 20:
        atoms = random.sample(mol.atoms, 20)
    else:
        atoms = mol.atoms
    atom_indices = [atom.index for atom in atoms]

    anagrad = -mol.calculate_forces()[atom_indices]
    numgrad = helpers.num_grad(mol,
                               mol.calculate_potential_energy,
                               atoms=atoms,
                               step=0.005*u.angstrom)
    assert (anagrad-numgrad).norm()/(3.0*len(atoms)) <= 5.0e-4 * u.eV / u.angstrom


@pytest.mark.parametrize('objkey', registered_types['hasmodel'])
def test_minimization_reduces_energy(objkey, request):
    mol = request.getfuncargvalue(objkey)
    e1 = mol.calculate_potential_energy()
    mol = request.getfuncargvalue(objkey)
    traj = mol.minimize()
    assert mol.calculate_potential_energy() < e1



@typedfixture('hasmodel')
def protein_stripped_lammps_amber_forcefield():
    mol = mdt.from_pdb('1yu8')
    mol = mdt.Molecule([res for res in mol.residues if res.type == "protein"])
    newmol = mdt.assign_forcefield(mol)
    newmol.set_energy_model(mdt.models.LAMMPSPotential)
    return newmol

@typedfixture('hasmodel')
def protein_lammps_amber_forcefield():
    mol = mdt.from_pdb('1yu8')
    newmol = mdt.assign_forcefield(mol)
    newmol.set_energy_model(mdt.models.LAMMPSPotential)
    return newmol

# TODO
def test_lammps_shake_constrain_hbonds(protein_stripped_lammps_amber_forcefield):
    mdmol = protein_stripped_lammps_amber_forcefield
    hydrogens = [atom for atom in mdmol.atoms if atom.atnum == 1]
    dist_before = []
    for atom in hydrogens:
        assert atom.num_bonds ==1
        d = atom.distance(atom.bonded_atoms[0])
        dist_before.append(d)
   
    mdmol.set_integrator(mdt.integrators.LAMMPSNvt, timestep=2.0*u.fs, temperature=300*u.kelvin, frame_interval=10.0*u.fs)
    mdmol.run(100.0*u.fs)
    
    numNonConstrained = 0
    for idx, atom in enumerate(hydrogens):
        assert atom.num_bonds ==1
        d = atom.distance(atom.bonded_atoms[0])
        if abs(d-dist_before[idx]).value_in(u.angstrom) > 0.0001:
            numNonConstrained += 1

    assert numNonConstrained < len(hydrogens)*0.1

# # TODO:
# def test_lammps_shake_constrain_water(protein_lammps_amber_forcefield):

