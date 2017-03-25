import moldesign as mdt

import pytest


def _assert_atom_indices_correct(m):
    for i, atom in enumerate(m.atoms):
        assert atom.index == i


def test_add_atom():
    atom = mdt.Atom('H')
    a2 = mdt.Atom('He')
    a3 = mdt.Atom('Li')
    mol = mdt.Molecule([atom])

    mol.atoms.append(a2)
    assert mol.num_atoms == 2
    _assert_atom_indices_correct(mol)

    a3.molecule = mol
    assert mol.num_atoms == 3
    _assert_atom_indices_correct(mol)


def test_add_atom_exceptions():
    atom = mdt.Atom(1)
    a2 = mdt.Atom(1)
    mol = mdt.Molecule([atom])
    m2 = mdt.Molecule([a2])

    with pytest.raises(ValueError):
        mol.atoms.append(atom)  # can't add same atom twice

    with pytest.raises(ValueError):
        mol.atoms.append(m2.atoms[0])  # can't add atom from another molecule

    assert mol.num_atoms == 1
    newatom = mdt.Atom(1)
    newatom.residue = mdt.Residue(name='TT')

    with pytest.raises(ValueError):
        mol.atoms.append(newatom)  # can't add atom that's part of a residue