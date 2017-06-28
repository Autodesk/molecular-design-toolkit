import pytest

import moldesign as mdt
from moldesign import units as u

from .test_ambertools_xface import protein_default_amber_forcefield, pdb1yu8
from .molecule_fixtures import *
from .test_qm_xfaces import h2_rhfwfn


def test_forcefield_atom_term_access(protein_default_amber_forcefield):
    mol = protein_default_amber_forcefield
    for atom in mol.atoms:
        assert atom.ff.ljepsilon.dimensionality == u.kcalpermol.dimensionality
        assert atom.ff.ljsigma.dimensionality == u.angstrom.dimensionality
        assert atom.ff.partial_charge.dimensionality == u.q_e.dimensionality

        assert atom.ff.pmdobj.idx == atom.index
        assert atom.ff.pmdobj.element == atom.atnum


@pytest.mark.screening
def test_forcefield_bond_term_access(protein_default_amber_forcefield):
    mol = protein_default_amber_forcefield
    for bond in mol.bonds:
        assert bond.ff.equilibrium_length.dimensionality == u.angstrom.dimensionality
        assert bond.ff.force_constant.dimensionality == (u.kcalpermol/u.angstrom).dimensionality

        pmdobj = bond.ff.pmdobj
        pmdatoms = sorted((pmdobj.atom1, pmdobj.atom2), key=lambda x:x.idx)
        assert (pmdatoms[0].idx, pmdatoms[1].idx) == (bond.a1.index, bond.a2.index)
        assert (pmdatoms[0].element, pmdatoms[1].element) == (bond.a1.atnum, bond.a2.atnum)


def test_atom_basis_function_returns_none_if_no_wfn(h2):
    for atom in h2.atoms:
        assert atom.basis_functions is None


def test_atom_ffterms_returns_none_if_no_ff(h2):
    for atom in h2.atoms:
        assert atom.ff is None


def test_bond_ffterms_returns_none_if_no_ff(h2):
    for bond in h2.bonds:
        assert bond.ff is None


def test_basis_function_atom_access(h2_rhfwfn):
    mol = h2_rhfwfn
    for atom in mol.atoms:
        assert len(atom.basis_functions) == 1  # good ol' sto-3g
        assert len(atom.basis_functions[0].primitives) == 3


@pytest.mark.screening
def test_atom_property_access_to_mulliken_charges(h2_rhfwfn):
    mol = h2_rhfwfn
    for atom in mol.atoms:
        assert abs(atom.properties.mulliken) <= 1e-5 * u.q_e


def test_atom_properties_are_empty_dict_if_nothings_computed(h2):
    empty = mdt.utils.DotDict()
    for atom in h2.atoms:
        assert atom.properties == empty
