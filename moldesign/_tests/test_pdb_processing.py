"""
This file collects tests relating to the many special cases we deal with when
trying to process biomolecular structures from the pdb.

These currently just test the specific cases that we've implemented fixes for
"""

import moldesign as mdt
from .helpers import get_data_path


def test_missing_terminal_atoms_3ac2():
    """ Tests that we can still detect terminal residues even if elements of the peptide bonds
    are missing from the structure
    """
    mol = mdt.read(get_data_path('3ac2.pdb'))
    assert mol.chains.A.n_terminal is not None
    assert mol.chains.A.c_terminal is not None


def test_disulfide_bond_detection_1hpk():
    mol = mdt.read(get_data_path('1hpk.pdb'))
    mol = mdt.interfaces.ambertools._prep_for_tleap(mol)
    for residx in (0, 21, 49, 61, 73, 78):
        residue = mol.residues[residx]
        assert residue.resname == 'CYX'
        assert residue.name[:3] == 'CYX'


def test_negative_residue_numbers_2jaj():
    mol = mdt.read(get_data_path('2jaj.pdb'))
    res = mol.chains.B.residues[0]
    assert res.pdbindex == -4
    assert res.index == 272
    assert res.name == 'GLY-4'

    mread = mdt.read(mol.write('pdb'),
                     format='pdb')
    res = mread.chains.B.residues[0]
    assert res.pdbindex == -4
    assert res.index == 272
    assert res.name == 'GLY-4'


def test_numeric_residue_name_1PYN():
    """ The ligand in this residue is named "941", which causes a little trickiness
    """
    import parmed

    mol = mdt.read(get_data_path('1pyn.pdb'))
    ligand = mol.residues[283]

    params = mdt.parameterize(ligand, charges='gasteiger')
    params.lib.put('/tmp/tmp.lib')

    contents = parmed.load_file('/tmp/tmp.lib')
    assert len(contents) == 1
    assert contents.keys()[0] == '941'





