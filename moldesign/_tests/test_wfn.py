import pytest
import numpy as np

import moldesign as mdt
from moldesign import units as u
from moldesign.interfaces.pyscf_interface import basis_values

from .molecule_fixtures import *
from . import helpers


__PYTEST_MARK__ = 'wfn'

TESTSYSTEMS = ['h2_rhf_augccpvdz', 'h2_rhf_sto3g', 'acetylene_dft_631g']


@pytest.mark.parametrize('molkey', TESTSYSTEMS)
def test_pyscf_orbital_grid_works(molkey, request):
    """ Tests the basic input/output of the pyscf basis_values function

    Doesn't actually test the values directly - just that the answers are mathematically consistent
    """
    mol = request.getfixturevalue(molkey)
    wfn = mol.wfn
    nbasis = len(wfn.aobasis)

    coords = u.array([mol.com,
                      np.zeros(3)*u.angstrom,
                      10.0 * np.ones(3) * u.angstrom,
                      np.ones(3)*u.nm])

    # First - check that the shape is appropriate if called without orbital coefficients
    values_nocoeffs = basis_values(mol, wfn.aobasis, coords)
    assert values_nocoeffs.shape == (len(coords), nbasis)
    assert (values_nocoeffs[-1] == values_nocoeffs[-2]).all()  # these 2 coordinates are the same

    # Second - explicitly send orbital coefficients for first 2 basis functions
    coeffs = np.zeros((2, nbasis))
    coeffs[:2, :2] = np.identity(2)
    vals_b0 = basis_values(mol, wfn.aobasis, coords, coeffs=coeffs)
    assert vals_b0.shape == (len(coords), len(coeffs))
    np.testing.assert_allclose(values_nocoeffs[:,:2], vals_b0)

    # Third - send symmetric and anti-symmetric combinations of basis functions and check answers
    plusminus = np.zeros((2, nbasis))
    plusminus[:2, :2] = 1.0 / np.sqrt(2)
    plusminus[1,1] = -1.0 / np.sqrt(2)
    vals_plusminus = basis_values(mol, wfn.aobasis, coords, coeffs=plusminus)
    assert vals_plusminus.shape == (len(coords), len(coeffs))
    helpers.assert_almost_equal(vals_plusminus[:,0],
                                (vals_b0[:,0] + vals_b0[:,1])/np.sqrt(2))
    helpers.assert_almost_equal(vals_plusminus[:,1],
                                (vals_b0[:,0] - vals_b0[:,1])/np.sqrt(2))


@pytest.mark.parametrize('molkey', TESTSYSTEMS)
def test_basis_function_3d_grids_same_in_pyscf_and_mdt(molkey, request):
    mol = request.getfixturevalue(molkey)
    randocoords = 6.0 * u.angstrom * (np.random.rand(200, 3) - 0.5)
    pyscf_vals = basis_values(mol, mol.wfn.aobasis, randocoords)
    with np.errstate(under='ignore'):
        mdt_vals = mol.wfn.aobasis(randocoords)
    helpers.assert_almost_equal(mdt_vals, pyscf_vals, decimal=6)


@pytest.mark.parametrize('molkey', ['h2_rhf_augccpvdz', 'h2_rhf_sto3g'])
@pytest.mark.screening
def test_pyscf_basis_function_space_integral_normalized(molkey, request):
    mol = request.getfixturevalue(molkey)
    grid = mdt.mathutils.padded_grid(mol.positions, 8.0 * u.angstrom, npoints=150)
    points = grid.allpoints()

    pyscf_vals = basis_values(mol, mol.wfn.aobasis, points)
    assert pyscf_vals.shape == (len(points), len(mol.wfn.aobasis))

    pyscf_normsqr = (pyscf_vals ** 2).sum(axis=0) * grid.dx * grid.dy * grid.dz
    helpers.assert_almost_equal(pyscf_normsqr, 1.0,
                                decimal=4)
