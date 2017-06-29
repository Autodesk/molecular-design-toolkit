import pytest

import moldesign as mdt
from .molecule_fixtures import *
from . import helpers

__PYTEST_MARK__ = 'internal'


def test_create_atom_with_element_as_name():
    he_plus = mdt.Atom("He", position=np.ones(3) * u.nm, formal_charge=1*u.q_e)
    assert he_plus.atomic_number == 2
    helpers.assert_almost_equal(he_plus.position, 10 * np.ones(3) * u.angstrom)
    assert he_plus.name == 'He'
    assert he_plus.element == 'He'
    assert he_plus.formal_charge == 1 * u.q_e


def test_create_atom_with_uppercase_element():
    cl_minus = mdt.Atom("CL2", element='CL', formal_charge=-1)
    assert cl_minus.formal_charge == -1 * u.q_e
    assert cl_minus.atnum == 17
    assert cl_minus.name == 'CL2'
    assert cl_minus.element == 'Cl'


def test_create_atom_with_unrelated_name():
    carbon_named_bob = mdt.Atom("bob", element='c')
    assert carbon_named_bob.name == 'bob'
    assert carbon_named_bob.atnum == 6
    assert carbon_named_bob.element == 'C'


def test_create_atom_with_element_name_as_first_char():
    carbon_alpha = mdt.Atom('ca')
    assert carbon_alpha.name == 'ca'
    assert carbon_alpha.atnum == 6
    assert carbon_alpha.element == 'C'


def test_create_atom_with_atomic_number():
    h = mdt.Atom('ca', atnum=1)
    assert h.atnum == 1
    assert h.element == 'H'
    assert h.name == 'ca'

    with pytest.raises(AssertionError):
        mdt.Atom('ca', atnum=1, element='He')


def test_add_atom(h2):
    newatom = mdt.Atom('Xe')
    h2.add_atom(newatom)
    assert newatom in h2.atoms
    assert h2.num_atoms == 3
    assert h2.num_residues == 1
    assert h2.num_chains == 1
    assert newatom.residue is h2.residues[0]
    assert newatom.chain is h2.chains[0]


def test_add_already_owned_atoms(h2):
    h2cpy = h2.copy()
    h2cpy.translate([10, 11, 12] * u.angstrom)

    h2.add_atoms(h2cpy.atoms)

    assert h2.num_atoms == 4
    assert h2cpy.num_atoms == 2

    for atom in h2cpy.atoms:
        assert atom.molecule is h2cpy
        assert atom.residue is h2cpy.residues[0]
        assert atom.chain is h2cpy.chains[0]

    for newatom in h2.atoms[2:]:
        assert newatom.molecule is h2
        assert newatom.residue is h2.residues[1]
        assert newatom.chain is h2.chains[1]

    helpers.assert_almost_equal(h2.positions[2:],
                                h2cpy.positions)


def test_atom_shortstr_3aid(pdb3aid):
    assert pdb3aid.atoms[10]._shortstr() == 'CA #10 in A.GLN2'


def test_atom_shortstr_benzene(benzene):
    for atom in benzene.atoms:
        assert atom._shortstr() == '%s #%s' % (atom.name, atom.index)


def test_atom_bond_iterators(benzene):
    atom = benzene.atoms[0]
    assert atom.atnum == 6  # sanity check
    atombonds = benzene.bond_graph[atom]

    assert set(atom.bonded_atoms) == set(atombonds.keys())
    assert atom.bond_graph == atombonds

    assert len(atom.heavy_bonds) == 2
    heavybonds = sorted(atom.heavy_bonds, key=lambda x:x.order)
    assert len(heavybonds) == 2
    assert heavybonds[0].order == 1 and heavybonds[1].order == 2

    assert len(atom.bonds) == 3
    hbond, csingle, cdouble = False, False, False
    for bond in atom.bonds:
        assert bond.a1 is atom  # since "atom" has the lowest index, it's always a1
        assert bond.partner(atom) is bond.a2

        if bond.a2.atnum == 1:
            assert hbond is False
            hbond = True
            assert bond.order == 1
            assert len(bond.a2.heavy_bonds) == 0
        else:
            assert bond.a2.atnum == 6
            if bond.order == 1:
                assert csingle is False
                csingle = True
                assert bond == heavybonds[0]
            else:
                assert cdouble is False
                cdouble = True
                assert bond.order == 2
                assert bond == heavybonds[1]

    assert hbond and csingle and cdouble





#def test_add_atom_to_residues(pdb3aid):
#    res = pdb3aid.residues[5]
#    newatom = mdt.Atom('Ta', residue=res)
#    pdb3aid.add_atom(newatom)

#    assert newatom.molecule is pdb3aid
#    assert newatom.chain is res.chain
#    assert newatom.residue is res

#    assert newatom in res