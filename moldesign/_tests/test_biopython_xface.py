
import pytest
import Bio.PDB
import moldesign as mdt

from .helpers import get_data_path
from .molecule_fixtures import pdb3aid


@pytest.fixture
def biopy_3aid():
    parser = Bio.PDB.PDBParser()
    structure = parser.get_structure('3aid', get_data_path('3aid.pdb'))
    return structure


def test_biopy_to_mdt(biopy_3aid, pdb3aid):
    mol = mdt.interfaces.biopython_to_mol(biopy_3aid)
    assert mol.num_atoms == pdb3aid.num_atoms
    assert mol.num_residues == pdb3aid.num_residues
    assert mol.numchains == pdb3aid.num_chains
