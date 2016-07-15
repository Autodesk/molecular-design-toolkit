import pytest


def test_openmm_retains_atom_numbering(fixture_name, request):
    mol = request.getfuncargvalue(fixture_name)
    assert mol.num_atoms == mol.energy_model.sim.system.getNumParticles()
    for iatom, (atom, stkatom) in enumerate(zip(mol.atoms,
                                                mol.energy_model.sim.topology.atoms())):
        assert atom.elem == stkatom.element.symbol
        assert iatom == atom.index == stkatom.index


def test_only_hydrogen_bonds_constrained(fixture_name, request):
    mol = request.getfuncargvalue(fixture_name)
    sim = mol.energy_model.sim
    all_hydrogens = set([atom for atom in mol.atoms if atom.atnum == 1])
    for ibond in xrange(mol.num_atoms):
        i, j, val = sim.system.getConstraintParameters(ibond)
        ai = mol.atoms[i]
        aj = mol.atoms[j]
        if ai.elem == 'H':  # exactly one of the atoms should be hydrogen
            assert aj.elem != 'H'
            h = ai
        else:
            assert aj.elem == 'H'
            h = aj
        assert h in all_hydrogens
        all_hydrogens.pop(h)
    assert len(all_hydrogens) == 0


