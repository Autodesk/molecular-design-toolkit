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

registered_types = {key:val[:] for key,val in molecule_standards.items()}

__all__ = ['registered_types']


def typedfixture(*types, **kwargs):
    """This is a decorator that lets us associate fixtures with one or more arbitrary types.
    We'll later use this type to determine what tests to run on the result"""

    def fixture_wrapper(func):
        __all__.append(func.__name__)
        for t in types:
            registered_types.setdefault(t, []).append(func.__name__)
        return pytest.fixture(**kwargs)(func)

    return fixture_wrapper


######################################
# Various python objects
TESTDICT = collections.OrderedDict((('a', 'b'),
                                    ('c', 3),
                                    ('d', 'e'),
                                    ('a', 1),
                                    (3, 35)))


@typedfixture('pickleable')
def dotdict():
    dd = utils.DotDict(TESTDICT)
    return dd


# Some objects with units
@typedfixture('pickleable')
def list_of_units():
    return [1.0 * u.angstrom, 1.0 * u.nm, 1.0 * u.a0]


@typedfixture('pickleable', 'equality')
def simple_unit_array():
    return np.array([1.0, -2.0, 3.5]) * u.angstrom


@typedfixture('pickleable', 'equality')
def unit_number():
    return 391.23948 * u.ureg.kg * u.angstrom / u.alpha


######################################
# Atom objects
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
# Hydrogen-related objects
@typedfixture('molecule')
def h2_harmonic(h2):
    mol = h2
    SPRING_CONSTANT = 1.0 * u.kcalpermol / (u.angstrom ** 2)
    model = mdt.models.HarmonicOscillator(k=SPRING_CONSTANT)
    integrator = mdt.integrators.VelocityVerlet(timestep=0.5*u.fs, frame_interval=30)
    mol.set_energy_model(model)
    mol.set_integrator(integrator)
    return mol


@typedfixture('pickleable')
def atom_bond_graph(h2):
    return h2.bond_graph[h2.atoms[0]]


@typedfixture('pickleable')
def mol_bond_graph(h2):
    return h2.bond_graph


@typedfixture('pickleable')
def mol_wfn(h2_rhf_sto3g):
    return h2_rhf_sto3g.copy().wfn


@typedfixture('pickleable')
def mol_properties(h2_rhf_sto3g):
    return h2_rhf_sto3g.copy().properties


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
pickleable = registered_types['pickleable'] + moldesign_objects

__all__.extend(['moldesign_objects', 'pickleable'])
