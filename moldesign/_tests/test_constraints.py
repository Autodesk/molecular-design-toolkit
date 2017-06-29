""" Tests constraint routines
"""
from builtins import range
import pytest
import numpy as np

import moldesign as mdt
from moldesign import units as u
from .molecule_fixtures import *
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
def four_particle_45_twist():
    mol = helpers._make_mol_with_n_hydrogens(4)
    mol.positions = u.nm*[[0.1, 0.0, -0.5],
                          [0.0, 0.0, -0.5],
                          [0.0, 0.0, 0.5],
                          [0.2, -0.2, 0.5]]

    for iatom in range(3):
        mol.atoms[iatom].bond_to(mol.atoms[iatom+1], 1)

    return mol


@typedfixture('constraint')
def dihedral_constraint_satisfied(four_particle_45_twist):
    mol = four_particle_45_twist
    c = mol.constrain_dihedral(*mol.atoms)
    return mol, c


@typedfixture('constraint')
def dihedral_constraint_unsatisfied(dihedral_constraint_satisfied):
    mol, c = dihedral_constraint_satisfied
    mol.atoms[0].z += 0.1*u.angstrom
    mol.atoms[1].x += -1.0 * u.angstrom
    mol.atoms[2].z -= -0.6 * u.angstrom
    return mol, c


@typedfixture('constraint')
def angle_constraint_satisfied(four_particle_45_twist):
    mol = four_particle_45_twist
    c = mol.constrain_angle(*mol.atoms[:3])
    return mol, c


@typedfixture('constraint')
def angle_constraint_unsatisfied(angle_constraint_satisfied):
    mol, c = angle_constraint_satisfied
    mol.atoms[1].x += -1.0 * u.angstrom
    mol.atoms[2].z -= -0.6 * u.angstrom
    return mol, c


@typedfixture('constraint')
def distance_constraint_satisfied(four_particle_45_twist):
    mol = four_particle_45_twist
    c = mol.constrain_distance(mol.atoms[1], mol.atoms[2])
    return mol, c


@typedfixture('constraint')
def distance_constraint_unsatisfied(distance_constraint_satisfied):
    mol, c = distance_constraint_satisfied
    mol.atoms[1].x += -1.0 * u.angstrom
    mol.atoms[2].z -= -0.6 * u.angstrom
    return mol, c


@typedfixture('constraint')
def atom_position_constraint_unsatisfied(four_particle_45_twist):
    mol = four_particle_45_twist
    c = mol.constrain_atom(mol.atoms[2])
    mol.atoms[2].x += 1.0*u.angstrom  # the gradient is singular when exactly satisfied
    return mol, c


@typedfixture('constraint')
def atom_coordinate_constraint_satisfied(four_particle_45_twist):
    mol = four_particle_45_twist
    c = mdt.geom.FixedCoordinate(mol.atoms[3], vector=np.array([1,1,-1]))
    return mol, c


@typedfixture('constraint')
def angle_constraint_unsatisfied(atom_coordinate_constraint_satisfied):
    mol, c = atom_coordinate_constraint_satisfied
    mol.atoms[3].x += -1.0 * u.angstrom
    mol.atoms[3].y -= -1.0 * u.angstrom
    return mol, c


@pytest.mark.screening
def test_distance_constraint(distance_constraint_satisfied):
    mol, c = distance_constraint_satisfied

    # satisfied
    np.testing.assert_allclose(c.current().value_in(u.nm),
                               1.0)
    assert c.satisfied()
    assert abs(c.error().value_in(u.nm)) <= 1e-10

    # unsatisfied
    mol.atoms[1].z = -0.6 * u.nm
    assert not c.satisfied()
    assert abs(c.error().value_in(u.nm) - 0.1) <= 1e-10


@pytest.mark.parametrize('objkey', registered_types['constraint'])
def test_constraint_gradient(objkey, request):
    mol, c = request.getfixturevalue(objkey)

    calc_grad = c.gradient()
    num_grad = helpers.num_grad(mol, c.error)
    np.testing.assert_allclose(num_grad.defunits_value(),
                               calc_grad.defunits_value(),
                               atol=5.0*helpers.DEFSTEP.defunits_value())


@pytest.mark.screening
def test_dihedral_constraint_errors(four_particle_45_twist):
    mol = four_particle_45_twist
    constraint = mol.constrain_dihedral(*mol.atoms)

    assert constraint.error() == 0.0
    assert constraint.satisfied()

    constraint.value = 30.0 * u.degrees
    np.testing.assert_allclose(constraint.error().value_in(u.degrees), 15.0)
    assert not constraint.satisfied()


    constraint.value = 60.0 * u.degrees
    np.testing.assert_allclose(constraint.error().value_in(u.degrees), -15.0)
    assert not constraint.satisfied()


def test_dihedral_constraint_errors_at_0(four_particle_45_twist):
    mol = four_particle_45_twist
    mol.atoms[3].position = [0.1, 0.0, 0.5] * u.angstrom

    constraint = mol.constrain_dihedral(*mol.atoms)
    assert constraint.error() == 0.0
    assert constraint.satisfied()

    for angle in (0, -10, 10, 90, -90) * u.degrees:
        constraint.value = angle
        np.testing.assert_allclose(constraint.error().value_in(u.degrees), -angle.magnitude)


def test_cannot_add_duplicate_constraints(mol_with_zerocharge_params):
    mol = mol_with_zerocharge_params
    mol.constrain_distance(*mol.atoms[:2])
    mol.constrain_angle(*mol.atoms[:3])
    mol.constrain_dihedral(*mol.atoms[:4])
    mol.constrain_atom(mol.atoms[0])
    mol.constrain_hbonds()

    # Failure on adding duplicates
    with pytest.raises(KeyError):
        mol.constrain_distance(*mol.atoms[:2])
    with pytest.raises(KeyError):
        mol.constrain_angle(*mol.atoms[:3])
    with pytest.raises(KeyError):
        mol.constrain_dihedral(*mol.atoms[:4])
    with pytest.raises(KeyError):
        mol.constrain_atom(mol.atoms[0])
    with pytest.raises(KeyError):
        mol.constrain_hbonds()

    # These should be ok
    mol.constrain_distance(*mol.atoms[1:3])
    mol.constrain_angle(*mol.atoms[1:4])
    mol.constrain_dihedral(*mol.atoms[1:5])
    mol.constrain_atom(mol.atoms[1])


def test_degrees_of_freedom_with_constraints(benzene, h2):
    water = mdt.from_smiles('[OH2]')
    water.residues[0].resname = 'HOH'

    benzene.residues[0].resname = 'UNL'
    h2.residues[0].resname = 'H2'

    mol = benzene.combine(h2, water, water)

    assert mol.num_atoms == 20  # assert I can count
    assert mol.dynamic_dof == 3 * mol.num_atoms

    full_dof = mol.dynamic_dof

    def _constrain_direction():
        c = mdt.geom.constraints.FixedCoordinate(mol.atoms[0], np.array([0.0, 0, 1.0]))
        mol.constraints.append(c)

    for constrainer, num_dof, args in ((mol.constrain_distance, 1, mol.atoms[:2]),
                                       (mol.constrain_angle, 1, mol.atoms[:3]),
                                       (mol.constrain_dihedral, 1, mol.atoms[:4]),
                                       (mol.constrain_atom, 3, mol.atoms[:1]),
                                       (mol.constrain_hbonds, 11, [True]),
                                       (_constrain_direction, 1, [])):
        constrainer(*args)
        assert mol.dynamic_dof == full_dof-num_dof
        mol.clear_constraints()
        assert mol.dynamic_dof == full_dof


def test_degrees_of_freedom_with_integrator_constraints(benzene):
    water = mdt.from_smiles('[OH2]')
    water.residues[0].resname = 'HOH'
    mol = benzene.combine(water)

    full_dof = mol.dynamic_dof
    assert full_dof == 3 * mol.num_atoms

    mol.set_energy_model(mdt.models.OpenMMPotential)
    mol.set_integrator(mdt.integrators.OpenMMLangevin)

    mol.integrator.params.remove_translation = False
    mol.integrator.params.remove_rotation = False
    mol.integrator.params.constrain_hbonds = False
    mol.integrator.params.constrain_water = False
    assert mol.dynamic_dof == full_dof

    mol.integrator.params.constrain_hbonds = True
    mol.integrator.params.constrain_water = False
    assert mol.dynamic_dof == full_dof - 8

    mol.integrator.params.constrain_hbonds = False
    mol.integrator.params.constrain_water = True
    assert mol.dynamic_dof == full_dof - 3

    mol.integrator.params.constrain_hbonds = True
    mol.integrator.params.constrain_water = True
    assert mol.dynamic_dof == full_dof - 9

    mol.clear_constraints()  # does nothing, constraints are in the integrator
    assert mol.dynamic_dof == full_dof - 9

    mol.integrator.params.constrain_hbonds = False
    mol.integrator.params.constrain_water = False
    assert mol.dynamic_dof == full_dof

    mol.integrator.params.remove_translation = True
    mol.integrator.params.remove_rotation = True
    mol.integrator.params.constrain_hbonds = True
    mol.integrator.params.constrain_water = True
    assert mol.dynamic_dof == full_dof - 9 - 6



def test_degrees_of_freedom_with_integrator_reorientation(pdb3aid):
    mol = pdb3aid

    full_dof = mol.dynamic_dof

    mol.set_energy_model(mdt.models.OpenMMPotential)
    mol.set_integrator(mdt.integrators.OpenMMLangevin)

    mol.integrator.params.remove_translation = False
    mol.integrator.params.remove_rotation = False
    mol.integrator.params.constrain_hbonds = False
    mol.integrator.params.constrain_water = False
    assert mol.dynamic_dof == full_dof

    mol.integrator.params.remove_translation = False
    mol.integrator.params.remove_rotation = True
    assert mol.dynamic_dof == full_dof - 3

    mol.integrator.params.remove_translation = True
    mol.integrator.params.remove_rotation = False
    assert mol.dynamic_dof == full_dof - 3

    mol.integrator.params.remove_translation = True
    mol.integrator.params.remove_rotation = True
    assert mol.dynamic_dof == full_dof - 6
