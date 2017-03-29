import pytest
import moldesign as mdt


@pytest.fixture(scope='function')
def hch():
    a1 = mdt.Atom(name='H1', atnum=1)
    a2 = mdt.Atom('C')
    a3 = mdt.Atom(name='H2', atnum=1)
    mol = mdt.Molecule([a1, a2, a3])
    mol.bonds.create(a1, a2, 1)
    a2.bonds.create(a3, 1)
    return mol

@pytest.fixture(scope='function')
def h2_atoms():
    a1 = mdt.Atom(name='H1', atnum=1)
    a2 = mdt.Atom(name='H2', atnum=1)
    return a1, a2


@pytest.fixture(scope='function')
def h2(h2_atoms):
    a1, a2 = h2_atoms
    a1.bonds.create(a2, 3)
    mol = mdt.Molecule([a1, a2])
    return mol


def test_bonds_outside_molecule(h2_atoms):
    a1, a2 = h2_atoms
    a1.bonds.create(a2, 3)
    assert a1 in a2.bonds
    assert a2 in a1.bonds


def test_create_molecule_with_prebonded_atoms(h2):
    mol = h2
    a1, a2 = mol.atoms
    assert mol.num_bonds == 1
    bonds = list(mol.bonds)
    assert len(bonds) == 1
    b = bonds[0]
    assert b.a1 is a1
    assert b.a2 is a2
    assert b.order == 3
    assert a1 in a2.bonds
    assert a2 in a1.bonds
    assert a1.bonds[a2].order == 3
    assert a1.bonds[a2].order == 3
    assert list(a1.bonds.iteritems()) == [(a2, 3)]
    assert list(a2.bonds.iteritems()) == [(a1, 3)]


def test_create_bonds_in_molecule(hch):
    mol = hch
    a1, a2, a3 = mol.atoms
    assert mol.num_bonds == 2
    assert a1 in a2.bonds
    assert a1 not in a3.bonds
    assert a3 not in a1.bonds

    assert a2 not in a2.bonds
    assert a1 not in a1.bonds
    assert a1 not in a3.bonds


def test_cant_bond_same_two_atoms_twice(hch):
    mol = hch
    a1, a2, a3 = mol.atoms

    with pytest.raises(ValueError):
        a3.bonds.create(a2, 3)

    with pytest.raises(ValueError):
        mol.bonds.create(a1, a2, 5)


def test_modify_bond_order(hch):
    mol = hch
    a1, a2, a3 = mol.atoms

    mol.bonds[a1, a2].order = 55
    assert a1.bonds[a2].order == 55
    for atom, order in a1.bonds.iteritems():
        if atom is a2:
            assert order == 55


def test_delete_bond(hch):
    mol = hch
    assert mol.num_bonds == 2
    bond = mol.atoms[0].bonds[mol.atoms[1]]
    assert bond in mol.bonds
    bond.delete()
    assert mol.num_bonds == 1
    assert mol.atoms[1] not in mol.atoms[0].bonds


def test_deleting_atom_deletes_its_bonds_from_the_molecule_only(hch):
    mol = hch
    atom = mol.atoms[0]
    nbr = mol.atoms[1]
    assert atom.num_bonds == 1
    assert nbr.num_bonds == 2
    assert mol.num_bonds == 2

    atom.molecule = None
    assert atom.num_bonds == 1  # it keeps its bonds, but they are not in the molecule any more
    assert mol.num_bonds == 1
    assert nbr.num_bonds == 1

    assert len(mol.bonds._graph) == 2
    for v in mol.bonds._graph.itervalues():
        assert len(v) == 1
