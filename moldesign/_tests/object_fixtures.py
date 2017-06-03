"""
This module is mostly _fixtures_ - functions that create various objects for downstream testing.

There are a few tests here too, mostly to make sure the objects are created correctly in the first
place.
"""

import collections

import numpy as np
import pytest

import moldesign as mdt
from .. import units as u
from .. import utils

from .molecule_fixtures import *

registered_types = {}
registered_types.update(molecule_standards)

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

TESTDICT = collections.OrderedDict((('a', 'b'),
                                    ('c', 3),
                                    ('d', 'e'),
                                    ('a', 1),
                                    (3, 35)))


@typedfixture('object')
def dotdict():
    dd = utils.DotDict(TESTDICT)
    return dd


# Some objects with units
@typedfixture('object')
def list_of_units():
    return [1.0 * u.angstrom, 1.0 * u.nm, 1.0 * u.a0]


@typedfixture('object', 'equality')
def simple_unit_array():
    return np.array([1.0, -2.0, 3.5]) * u.angstrom


@typedfixture('object', 'equality')
def unit_number():
    return 391.23948 * u.ureg.kg * u.ang / u.alpha


######################################
# Test underlying elements
@typedfixture('atom')
def carbon_atom():
    atom1 = mdt.Atom('C')
    return atom1


def test_carbon_atom(carbon_atom):
    assert carbon_atom.symbol == 'C'
    assert carbon_atom.mass == 12.0 * u.amu


@typedfixture('atom')
def carbon_copy(carbon_atom):
    atoms = carbon_atom.copy()
    return atoms


######################################
# Tests around hydrogen


@typedfixture('molecule')
def h2_harmonic(h2):
    mol = h2
    SPRING_CONSTANT = 1.0 * u.kcalpermol / (u.angstrom ** 2)
    model = mdt.models.HarmonicOscillator(k=SPRING_CONSTANT)
    integrator = mdt.integrators.VelocityVerlet(timestep=0.5*u.fs, frame_interval=30)
    mol.set_energy_model(model)
    mol.set_integrator(integrator)
    return mol


@typedfixture('trajectory')
def h2_trajectory(h2_harmonic):
    mol = h2_harmonic
    mol.atoms[0].x = 1.0 * u.angstrom
    mol.momenta *= 0.0
    traj = mol.run(500)
    return traj


@typedfixture('molecule')
def h2_traj_tempmol(h2_trajectory):
    return h2_trajectory._tempmol


@typedfixture('molecule')
def h2_harmonic_copy(h2_harmonic):
    return mdt.Molecule(h2_harmonic)


@typedfixture('submolecule')
def copy_atoms_from_h2_harmonic(h2_harmonic):
    atoms = h2_harmonic.atoms.copy()
    return atoms

@typedfixture('molecule')
def h2_harmonic_thats_been_copied(h2_harmonic):
    temp = mdt.Molecule(h2_harmonic)
    return h2_harmonic


@typedfixture('submolecule')
def h2_harmonic_atoms(h2_harmonic):
    return h2_harmonic.atoms


moldesign_objects = (registered_types['molecule'] +
                     registered_types['submolecule'] +
                     registered_types['trajectory'] +
                     registered_types['atom'])
all_objects = registered_types['object'] + moldesign_objects
