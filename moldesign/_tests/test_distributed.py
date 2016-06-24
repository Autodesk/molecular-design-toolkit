import moldesign as mdt
from moldesign import units as u


def test_amber_tleap():
    mol = mdt.read('../data/nuc.pdb')
    model = moldesign.methods.models.OpenMMPotential(implicit_solvent='obc',
                                                     cutoff=8.0*u.angstrom,
                                                     force_field='ff14SB')
    mol.set_energy_model(model)
    result = mol.calculate(wait=True)
    assert 'potential_energy' in result and 'forces' in result