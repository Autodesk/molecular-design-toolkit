import pytest
import moldesign as mdt

from .helpers import check_topology_consistency


@pytest.fixture(scope='function')
def methane():
    mol = mdt.from_smiles('C')
    assert mol.num_atoms == 5
    return mol


def _check_added(mol, atom):
    assert atom in mol
    assert atom in mol.atoms
    assert atom is mol.atoms[atom.index]
    assert atom.residue.molecule is mol
    assert atom.residue is mol.residues[atom.residue.index]
    assert atom.residue[atom.name] is atom


def _check_isolated(atom, mol):
    assert atom.molecule is atom.chain is atom.residue is None
    assert atom not in mol


def test_mol_init_from_residues():
    atoms = mdt.AtomList()
    resa = mdt.Residue(name='A')
    resb = mdt.Residue(name='B')
    for i in xrange(1, 6):
        a1 = mdt.Atom(atnum=i, name='A%d'%i)
        a2 = mdt.Atom(atnum=i, name='B%d'%i)
        a1.residue = resa
        resb.add(a2)
        atoms.extend([a1, a2])

    assert resa.num_atoms == 5
    assert resb.num_atoms == 5

    mol = mdt.Molecule(atoms)
    check_topology_consistency(mol)

    assert mol.num_atoms == 10
    assert mol.num_residues == 2
    assert mol.num_chains == 1
    assert mol.chains[0].num_residues == 2

    for a1, a2 in zip(atoms, mol.atoms):
        assert a1 is a2


def test_add_atom_with_append(methane):
    mol = methane
    a2 = mdt.Atom('He')
    mol.atoms.append(a2)
    assert mol.num_atoms == 6
    _check_added(mol, a2)
    check_topology_consistency(mol)
    assert mol.num_residues == 1
    assert mol.num_chains == 1


def test_add_atom_by_assigning_molecule(methane):
    mol = methane
    a3 = mdt.Atom('Li')
    a3.molecule = mol
    assert mol.num_atoms == 6
    _check_added(mol, a3)
    check_topology_consistency(mol)
    assert mol.num_residues == 1
    assert mol.num_chains == 1


def test_add_atom_in_another_residue(methane):
    mol = methane
    a4 = mdt.Atom('Ar')
    a4.residue = mdt.Residue(name='Be')
    mol.atoms.append(a4)
    _check_added(mol, a4)
    check_topology_consistency(mol)
    assert mol.num_atoms == 6
    assert mol.num_residues == 2
    assert mol.num_chains == 1

def test_remove_atom_with_moleculelist(methane):
    mol = methane
    hatom = mol.atoms[3]
    mol.atoms.remove(hatom)
    assert mol.num_atoms == 4
    _check_isolated(hatom, mol)
    check_topology_consistency(mol)


def test_remove_atom_with_null_residue(methane):
    mol = methane
    hatom = mol.atoms[3]
    hatom.residue = None
    assert mol.num_atoms == 4
    _check_isolated(hatom, mol)
    check_topology_consistency(mol)


def test_remove_atom_with_null_molecule(methane):
    mol = methane
    hatom = mol.atoms[3]
    hatom.molecule = None
    assert mol.num_atoms == 4
    _check_isolated(hatom, mol)
    check_topology_consistency(mol)


def test_change_atom_name(methane):
    mol = methane
    atom = mol.atoms[0]
    atom.name = 'kevin'
    assert mol.residues[0]['kevin'] is atom
    check_topology_consistency(mol)


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
    newatom.residue.add(mdt.Atom(3))

    with pytest.raises(ValueError):
        mol.atoms.append(newatom)  # can't add atom that's part of a residue


def test_move_atom_between_residues(methane):
    mol = methane.combine(methane.copy())
    check_topology_consistency(mol)
    assert mol.num_residues == 2
    assert mol.num_chains == 2

    assert mol.atoms[0].residue is mol.residues[0]
    assert mol.atoms[0] in mol.residues[0]

    with pytest.raises(ValueError):
        # can't do it this way because it already belongs to residue 0
        mol.residues[1].add(mol.atoms[0])

    moveatom = mol.atoms[0]
    moveatom.residue = mol.residues[1]
    check_topology_consistency(mol)

    assert mol.residues[0].num_atoms == 4
    assert mol.residues[1].num_atoms == 6
    assert moveatom not in mol.residues[0]
    assert moveatom in mol.residues[1]


def test_move_atom_between_molecules(methane):
    m = methane
    m2 = methane.copy()

    with pytest.raises(ValueError):
        m2.atoms.append(m.atoms[0])

    moveatom = m.atoms[0]
    for bond in list(moveatom.bonds):
        bond.delete()

    moveatom.molecule = m2
    check_topology_consistency(m)
    check_topology_consistency(m2)
    assert moveatom.num_bonds == 0
    assert moveatom.molecule is m2
    assert m2.num_atoms == 6
    assert m2.num_chains == 2
    assert m2.num_residues == 2
    assert moveatom.residue.num_atoms == 1
    assert moveatom.chain.num_atoms == moveatom.chain.num_residues == 1


def test_move_residue_between_molecules(methane):
    mol = methane
    m2 = methane.copy()
    m2.residues[0].molecule = mol
    check_topology_consistency(mol)
    check_topology_consistency(m2)

    assert m2.num_residues == 0
    assert m2.num_atoms == 0
    assert mol.num_residues == 2
    assert mol.num_chains == 1  # because the residue is assigned to the default chain in both cases
    assert mol.num_atoms == 10


def test_insert_residue_into_chain_in_another_molecule(methane):
    mol = methane.combine(methane.copy())
    m3 = methane.copy()
    check_topology_consistency(m3)

    assert mol.num_residues == mol.num_chains == 2

    m3.residues[0].chain = mol.chains[0]
    check_topology_consistency(mol)
    check_topology_consistency(m3)

    assert mol.chains[0].num_residues == 2
    assert mol.chains[0].num_atoms == 10
    assert mol.num_atoms == 15
    assert mol.num_residues == 3
    assert mol.num_chains == 2


def test_insert_residue_into_chain_in_same_molecule(methane):
    mol = methane.combine(methane.copy())

    assert mol.num_residues == mol.num_chains == 2

    mol.chains[1].residues[0].chain = mol.chains[0]
    check_topology_consistency(mol)

    assert mol.chains[0].num_residues == 2
    assert mol.chains[0].num_atoms == 10
    assert mol.num_atoms == 10
    assert mol.chains[1].num_residues == 0
    assert mol.chains[1].num_atoms == 0


def test_remove_residue_explicitly(methane):
    mol = methane.combine(methane.copy())
    mol.residues.remove(mol.residues[1])
    check_topology_consistency(mol)

    assert mol.num_atoms == 5
    assert mol.num_residues == 1
    assert mol.num_chains == 2


def test_remove_residue_by_unassigning_molecule(methane):
    mol = methane.combine(methane.copy())
    mol.residues[1].molecule = None
    check_topology_consistency(mol)

    assert mol.num_atoms == 5
    assert mol.num_residues == 1
    assert mol.num_chains == 2


def test_remove_residue_by_unassigning_chain(methane):
    mol = methane.combine(methane.copy())
    mol.residues[1].chain = None
    check_topology_consistency(mol)

    assert mol.num_atoms == 5
    assert mol.num_residues == 1
    assert mol.num_chains == 2


def test_remove_chain_explicitly(methane):
    mol = methane.combine(methane.copy())
    mol.chains.remove(mol.chains[0])
    check_topology_consistency(mol)

    assert mol.num_atoms == 5
    assert mol.num_chains == 1
    assert mol.num_residues == 1


def test_remove_chain_by_unassigning_molecule(methane):
    m = methane
    m2 = methane.copy()
    m2.chains[0].name = 'C'
    mol = m.combine(m2)
    assert mol.num_residues == 2
    check_topology_consistency(mol)

    mol.chains[0].molecule = None
    check_topology_consistency(mol)
    assert mol.num_chains == 1


def test_move_chain_between_molecules_by_assignment(methane):
    m = methane.combine(methane.copy())
    m2 = methane.copy()

    m.chains[1].molecule = m2
    check_topology_consistency(m)
    check_topology_consistency(m2)

    assert m.num_atoms == 5
    assert m2.num_atoms == 10
    assert m.num_chains == 1
    assert m2.num_chains == 2
    assert m.num_residues == 1
    assert m2.num_residues == 2
