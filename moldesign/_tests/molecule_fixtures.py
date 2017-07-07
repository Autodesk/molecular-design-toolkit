import random
import pytest
import numpy as np

import moldesign as mdt
import moldesign.units as u
from moldesign.utils import exports
from .helpers import get_data_path


__all__ = ['molecule_standards']

molecule_standards = {}

def typedfixture(*types, **kwargs):
    """This is a decorator that lets us associate fixtures with one or more arbitrary types.
    We'll later use this type to determine what tests to run on the result"""

    def fixture_wrapper(func):
        __all__.append(func.__name__)
        for t in types:
            molecule_standards.setdefault(t, []).append(func.__name__)
        return pytest.fixture(**kwargs)(func)

    return fixture_wrapper


######################################
# Tests around PDB ID 3AID
@typedfixture('molecule')
def pdb3aid():
    mol = mdt.read(get_data_path('3aid.pdb'))
    return mol


@typedfixture('submolecule')
def ligand_residue_3aid(pdb3aid):
    unknown = pdb3aid.chains['A'](type='unknown')
    assert len(unknown) == 1
    return unknown[0]


@typedfixture('submolecule')
def ligand_3aid_atoms(ligand_residue_3aid):
    return ligand_residue_3aid.atoms


@typedfixture('molecule')
def ligand3aid(ligand_residue_3aid):
    newmol = mdt.Molecule(ligand_residue_3aid)
    return newmol


@typedfixture('molecule')
def ethylene_waterbox_2na_2cl():
    mol = mdt.from_smiles('C=C')
    solvated = mdt.add_water(mol, padding=15.0*u.angstrom, ion_concentration=0.6*u.molar)
    return solvated


@exports
@pytest.fixture
def random_atoms_from_3aid(pdb3aid):
    atoms = mdt.molecules.atomcollections.AtomList(random.sample(pdb3aid.atoms, 10))
    return atoms


@exports
@pytest.fixture(scope='session')
def cached_small_molecule():
    mol = mdt.from_smiles('CNCOS(=O)C')
    mol.positions += 0.001*np.random.random(mol.positions.shape)*u.angstrom  # move out of minimum
    return mol


@exports
@pytest.fixture
def small_molecule(cached_small_molecule):
    return cached_small_molecule.copy()


@exports
@pytest.fixture(scope='session')
def cached_benzene():
    return mdt.from_smiles('c1ccccc1')


@exports
@pytest.fixture
def benzene(cached_benzene):
    return cached_benzene.copy()


@typedfixture('molecule')
def h2():
    mol = mdt.Molecule([mdt.Atom('H1'),
                        mdt.Atom('H2')])
    mol.atoms[0].x = 0.5 * u.angstrom
    mol.atoms[1].x = -0.25 * u.angstrom
    mol.atoms[0].bond_to(mol.atoms[1], 1)
    return mol


@typedfixture('molecule')
def heh_plus():
    mol = mdt.Molecule([mdt.Atom('H'),
                        mdt.Atom('He')])
    mol.atoms[1].z = 1.0 * u.angstrom
    mol.charge = 1 * u.q_e
    mol.atoms[0].bond_to(mol.atoms[1], 1)
    return mol


@exports
@pytest.fixture(scope='session')
def cached_ethylene():
    return mdt.from_smiles('C=C')


@exports
@pytest.fixture
def ethylene(cached_ethylene):
    return cached_ethylene.copy()


@exports
@pytest.fixture(scope='session')
def cached_pdb1yu8():
    return mdt.read(get_data_path('1yu8.pdb'))

@exports
@pytest.fixture
def pdb1yu8():
    return mdt.read(get_data_path('1yu8.pdb'))


@exports
@pytest.fixture(scope='session')
def cached_mol_from_xyz():
    return mdt.read("""43
    c1.pdb
    C         -1.21700        1.04300        2.45300
    C         -0.14200        0.18700        2.19500
    C         -0.31600       -0.99500        1.46200
    C         -1.59800       -1.33100        1.02200
    C         -2.68200       -0.48500        1.28100
    C         -2.50400        0.70500        1.98200
    O         -3.53000        1.57400        2.25000
    O         -1.13200        2.21500        3.14700
    C          0.14200        2.61500        3.63500
    C          0.86700       -1.90600        1.12900
    C          1.10600       -1.99700       -0.40500
    O          2.06900       -1.52700        1.78600
    O          1.81300       -0.81200       -0.83600
    C          1.98100       -3.18300       -0.81100
    O          2.18000       -3.27400       -2.21400
    C          1.18300        0.32500       -1.28900
    C          0.13100        0.33800       -2.23300
    C         -0.38800        1.56800       -2.65400
    C          0.15900        2.77500       -2.20100
    C          1.23000        2.75900       -1.31000
    C          1.72800        1.53500       -0.85900
    O         -0.31100       -0.88000       -2.70600
    C         -1.43500       -0.90100       -3.58200
    H          0.84400        0.42500        2.57000
    H         -1.77600       -2.25500        0.47700
    H         -3.67600       -0.75700        0.93200
    H         -4.35900        1.18400        1.93700
    H         -0.03000        3.54300        4.18500
    H          0.84500        2.80900        2.81200
    H          0.56300        1.86400        4.31600
    H          0.63600       -2.91500        1.49500
    H          0.14200       -2.06200       -0.91500
    H          2.51900       -0.89500        1.20400
    H          1.52900       -4.10700       -0.41600
    H          2.97500       -3.07900       -0.36800
    H          1.34400       -3.10100       -2.67400
    H         -1.20600        1.60000       -3.36300
    H         -0.25300        3.71600       -2.56000
    H          1.68000        3.68900       -0.97300
    H          2.55900        1.49100       -0.16100
    H         -1.60900       -1.95000       -3.81500
    H         -1.22800       -0.35800       -4.51600
    H         -2.32500       -0.48000       -3.09700
    """, format='xyz')


@exports
@pytest.fixture
def mol_from_xyz(cached_mol_from_xyz):
    return cached_mol_from_xyz.copy()


@exports
@pytest.fixture(scope='session')
def cached_mol_from_sdf():
    return mdt.read("""
 OpenBabel02271712493D

 11 10  0  0  0  0  0  0  0  0999 V2000
    0.9114   -0.0615    0.0032 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4168    1.3750    0.0264 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9114    2.1503   -1.1831 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5568   -0.5828   -0.8916 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.5460   -0.6043    0.8807 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0053   -0.0987    0.0111 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6785    1.3847    0.0446 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7625    1.8663    0.9425 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.5568    1.6932   -2.1123 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.5460    3.1815   -1.1499 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0053    2.1774   -1.2097 H   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
  3  2  1  0  0  0  0
  4  1  1  0  0  0  0
  5  1  1  0  0  0  0
  6  1  1  0  0  0  0
  7  2  1  0  0  0  0
  8  2  1  0  0  0  0
  9  3  1  0  0  0  0
 10  3  1  0  0  0  0
 11  3  1  0  0  0  0
M  END
$$$$
""", format='sdf')

@exports
@pytest.fixture
def mol_from_sdf(cached_mol_from_sdf):
    return cached_mol_from_sdf.copy()


@typedfixture('molecule')
def nucleic():
    # ACTG.pdb contains a molecule generated using mdt.build_dna('ACTG')
    mol = mdt.read(get_data_path('ACTG.pdb'))
    return mol


########################################################################################
# Molecules with forcefields assigned - these use a session-scoped constructor w/ a copy factory

@exports
@pytest.fixture(scope='session')
def cached_mol_parameterized_with_zeros(cached_small_molecule):
    return _param_small_mol(cached_small_molecule.copy(), 'zero')


@typedfixture('hasmodel')
def mol_with_zerocharge_params(cached_mol_parameterized_with_zeros):
    return cached_mol_parameterized_with_zeros.copy()


@exports
@pytest.fixture(scope='session')
def cached_mol_parameterized_with_am1bcc(cached_small_molecule):
    """ We don't use this fixture directly, rather use another fixture that copies these results
    so that we don't have to repeatedly call tleap/antechamber
    """
    return _param_small_mol(cached_small_molecule.copy(), 'am1-bcc')


@typedfixture('hasmodel')
def mol_with_am1bcc_params(cached_mol_parameterized_with_am1bcc):
    return cached_mol_parameterized_with_am1bcc.copy()


@exports
@pytest.fixture(scope='session')
def cached_mol_parameterized_with_gasteiger(cached_small_molecule):
    """ We don't use this fixture directly, rather use another fixture that copies these results
    so that we don't have to repeatedly call tleap/antechamber
    """
    mol = cached_small_molecule.copy()
    mol.set_energy_model(mdt.models.GaffSmallMolecule, partial_charges='gasteiger')
    mol.energy_model.prep()
    return mol


@typedfixture('hasmodel')
def mol_with_gast_params(cached_mol_parameterized_with_gasteiger):
    return cached_mol_parameterized_with_gasteiger.copy()


def _param_small_mol(cached_small_molecule, chargemodel):
    mol = cached_small_molecule.copy()
    params = mdt.create_ff_parameters(mol, charges=chargemodel, baseff='gaff2')
    params.assign(mol)
    mol.set_energy_model(mdt.models.ForceField)
    return mol


@exports
@pytest.fixture(scope='session')
def cached_protein_with_default_amber_ff(cached_pdb1yu8):
    """ We don't use this fixture directly, rather use another fixture that copies these results
    so that we don't have to repeatedly call tleap/antechamber
    """
    mol = cached_pdb1yu8
    ff = mdt.forcefields.DefaultAmber()
    newmol = ff.create_prepped_molecule(mol)
    newmol.set_energy_model(mdt.models.ForceField)
    return newmol


@typedfixture('hasmodel')
def protein_default_amber_forcefield(cached_protein_with_default_amber_ff):
    return cached_protein_with_default_amber_ff.copy()


@exports
@pytest.fixture(scope='session')
def cached_h2_rhf_sto3g():
    mol = h2()  # fixture is not cached, so just call it directly
    mol.set_energy_model(mdt.models.PySCFPotential, basis='sto-3g', theory='rhf')
    mol.calculate(requests=['forces'])
    return mol


@exports
@pytest.fixture
def h2_rhf_sto3g(cached_h2_rhf_sto3g):
    return cached_h2_rhf_sto3g.copy()


@exports
@pytest.fixture(scope='session')
def cached_h2_rhf_augccpvdz():
    mol = h2()
    mol.set_energy_model(mdt.models.RHF, basis='aug-cc-pvdz')
    mol.calculate()
    return mol


@exports
@pytest.fixture
def h2_rhf_augccpvdz(cached_h2_rhf_augccpvdz):
    return cached_h2_rhf_augccpvdz.copy()


@exports
@pytest.fixture(scope='session')
def cached_acetylene_dft_631g():
    mol = mdt.from_smiles('C#C')
    mol.set_energy_model(mdt.models.B3LYP, basis='6-31g')
    mol.calculate()
    return mol


@exports
@pytest.fixture
def acetylene_dft_631g(cached_acetylene_dft_631g):
    return cached_acetylene_dft_631g.copy()
