import moldesign as mdt

import pytest

def _test_topology_consistency(mol):
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
    _test_topology_consistency(mol)

    assert mol.num_atoms == 10
    assert mol.num_residues == 2
    assert mol.num_chains == 1
    assert mol.chains[0].num_residues == 2


def test_add_atom():
    atom = mdt.Atom('H')
    a2 = mdt.Atom('He')
    a3 = mdt.Atom('Li')
    mol = mdt.Molecule([atom])
    assert mol.num_atoms == 1
    assert mol.atoms[0] is atom

    mol.atoms.append(a2)
    assert mol.num_atoms == 2
    _test_topology_consistency(mol)

    a3.molecule = mol
    assert mol.num_atoms == 3
    _test_topology_consistency(mol)


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