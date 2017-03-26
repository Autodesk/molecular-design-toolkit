import moldesign as mdt

import pytest


@pytest.fixture(scope='function')
def methane():
    mol = mdt.from_smiles('C')
    assert mol.num_atoms == 5
    return mol


def _check_topology_is_consistent(mol):
    all_residues = set(mol.residues)
    all_chains = set(mol.chains)
    found_chains = set(atom.chain for atom in mol.atoms)
    found_residues = set(atom.residue for atom in mol.atoms)

    assert len(all_chains) == len(mol.chains)
    assert all_chains == found_chains
    assert len(found_residues) == len(all_residues)

    residues_from_chains = set()
    for ichain, chain in enumerate(mol.chains):
        assert chain.index == ichain
        assert mol.chains[chain.index] is chain
        if chain.name is not None:
            assert mol.chains[chain.name] is chain
        assert chain.molecule is mol

        for residue in chain:
            assert residue not in residues_from_chains, 'Residue in more than 1 chain'
            residues_from_chains.add(residue)
            assert residue in all_residues
            assert residue.molecule is mol
            assert residue.chain is chain
            assert mol.residues[residue.index] is residue
            if residue.name is not None:
                assert chain[residue.name] is residue

            for atom in residue:
                assert atom.residue is residue
                assert atom.chain is chain
                assert atom.molecule is mol
                assert mol.atoms[atom.index] is atom
                assert residue[atom.name] is atom

    assert residues_from_chains == all_residues


def _check_added(mol, atom):
    assert atom in mol
    assert atom in mol.atoms
    assert atom is mol.atoms[atom.index]
    assert atom.residue.molecule is mol
    assert atom.residue is mol.residues[atom.residue.index]
    assert atom.residue[atom.name] is atom


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
    _check_topology_is_consistent(mol)

    assert mol.num_atoms == 10
    assert mol.num_residues == 2
    assert mol.num_chains == 1
    assert mol.chains[0].num_residues == 2


def test_add_atom_with_append(methane):
    mol = methane
    a2 = mdt.Atom('He')
    mol.atoms.append(a2)
    assert mol.num_atoms == 6
    _check_added(mol, a2)
    _check_topology_is_consistent(mol)
    assert mol.num_residues == 1
    assert mol.num_chains == 1


def test_add_atom_by_assigning_molecule(methane):
    mol = methane
    a3 = mdt.Atom('Li')
    a3.molecule = mol
    assert mol.num_atoms == 6
    _check_added(mol, a3)
    _check_topology_is_consistent(mol)
    assert mol.num_residues == 1
    assert mol.num_chains == 1


def test_add_atom_in_another_residue(methane):
    mol = methane
    a4 = mdt.Atom('Ar')
    a4.residue = mdt.Residue(name='Be')
    mol.atoms.append(a4)
    _check_added(mol, a4)
    _check_topology_is_consistent(mol)
    assert mol.num_atoms == 6
    assert mol.num_residues == 2
    assert mol.num_chains == 1

def _check_isolated(atom, mol):
    assert atom.molecule is atom.chain is atom.residue is None
    assert atom not in mol


def test_remove_atom_with_moleculelist(methane):
    mol = methane
    hatom = mol.atoms[3]
    mol.atoms.remove(hatom)
    assert mol.num_atoms == 4
    _check_isolated(hatom, mol)
    _check_topology_is_consistent(mol)


def test_remove_atom_with_null_residue(methane):
    mol = methane
    hatom = mol.atoms[3]
    hatom.residue = None
    assert mol.num_atoms == 4
    _check_isolated(hatom, mol)
    _check_topology_is_consistent(mol)


def test_remove_atom_with_null_molecule(methane):
    mol = methane
    hatom = mol.atoms[3]
    hatom.molecule = None
    assert mol.num_atoms == 4
    _check_isolated(hatom, mol)
    _check_topology_is_consistent(mol)


def test_change_atom_name(methane):
    mol = methane
    atom = mol.atoms[0]
    atom.name = 'kevin'
    assert mol.residues[0]['kevin'] is atom
    _check_topology_is_consistent(mol)


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