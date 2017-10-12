import io
import pytest
import Bio.PDB
import moldesign as mdt

from .helpers import get_data_path
from .molecule_fixtures import pdb3aid


@pytest.fixture
def biopy_3aid():
    import gzip
    parser = Bio.PDB.PDBParser()
    s = gzip.open(get_data_path('3aid.pdb.gz'), 'r').read().decode('utf-8')
    structure = parser.get_structure('3aid', io.StringIO(s))
    return structure


@pytest.mark.screening
def test_biopy_to_mdt(biopy_3aid, pdb3aid):
    mol = mdt.interfaces.biopython_to_mol(biopy_3aid)
    assert mol.num_atoms == pdb3aid.num_atoms
    assert mol.num_residues == pdb3aid.num_residues
    assert mol.numchains == pdb3aid.num_chains
