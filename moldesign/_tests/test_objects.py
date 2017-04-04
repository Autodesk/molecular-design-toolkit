import itertools

import numpy as np
import pytest

from moldesign import units as u
from .object_fixtures import *


@pytest.mark.parametrize('objkey', 'dotdict ordered_dotdict'.split())
def test_dotdict(objkey, request):
    dd = request.getfuncargvalue(objkey)

    # basic lookups and functions
    assert len(dd) == len(TESTDICT)
    assert dd.a == TESTDICT['a']
    assert dd.d == TESTDICT['d']
    assert set(dd.keys()) == set(TESTDICT.keys())
    assert set(dd.values()) == set(TESTDICT.values())

    # __contains__ methods
    for item in TESTDICT:
        assert item in dd

    # item deletion
    del dd.d
    assert 'd' not in dd
    assert len(dd) == len(TESTDICT) - 1
    del dd[3]
    assert 3 not in dd

    # equivalence of items and attrs
    dd['k'] = 12345
    assert getattr(dd, 'k') == 12345
    setattr(dd, 'newkey', -42)
    assert dd['newkey'] == -42


def test_ordered_dotdict(ordered_dotdict):
    assert ordered_dotdict.keys() == TESTDICT.keys()
    assert ordered_dotdict.values() == TESTDICT.values()
    assert ordered_dotdict.items() == TESTDICT.items()

    for iterator_method in '__iter__ iterkeys itervalues iteritems'.split():
        odd_iter = getattr(ordered_dotdict, iterator_method)()
        test_iter = getattr(TESTDICT, iterator_method)()
        for gotval, testval in itertools.izip_longest(odd_iter, test_iter, fillvalue=None):
            assert gotval == testval, 'Iterator %s failed' % iterator_method


def test_h2_positions(h2):
    atom1, atom2 = h2.atoms
    assert (atom1.position == np.array([0.5, 0.0, 0.0]) * u.angstrom).all()
    assert atom2.x == -0.5 * u.angstrom
    assert atom1.distance(atom2) == 1.0 * u.angstrom


def test_h2(h2):
    mol = h2
    assert mol.num_atoms == 2
    assert mol.atoms[0].symbol == mol.atoms[1].symbol == 'H'


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