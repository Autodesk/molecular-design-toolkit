import pytest

import moldesign as mdt


@pytest.fixture(scope='function')
def mol_from_xyz():
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


@pytest.fixture(scope='function')
def mol_from_sdf():
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


@pytest.fixture(scope='function')
def mol_from_pdb():
    mol = mdt.from_pdb('3aid')
    return mdt.Molecule(mol.chains['A'].residues['ARQ401'])


@pytest.fixture(scope='function')
def mol_from_smiles():
    return mdt.from_smiles('CC')


@pytest.mark.parametrize('objkey',
                         'mol_from_smiles mol_from_pdb mol_from_xyz mol_from_sdf'.split())
def test_parameterization_from_formats(objkey, request):
    mol = request.getfuncargvalue(objkey)
    assert not mol.ff
    params = mdt.parameterize(mol, charges='gasteiger')
    _test_succesful_parameterization(mol, params)


@pytest.mark.parametrize('chargemodel',
                         'esp gasteiger zero am1-bcc'.split())
def test_charge_models(mol_from_smiles, chargemodel):
    mol = mol_from_smiles

    if chargemodel == 'esp':
        pytest.xfail("ESP not yet implemented")
    assert not mol.ff
    params = mdt.parameterize(mol, charges=chargemodel)
    _test_succesful_parameterization(mol, params)


def _test_succesful_parameterization(mol, params):
    assert params is not None
    assert mol.ff
    mol.set_energy_model(mdt.models.ForceField)
    mol.calculate()
    assert 'potential_energy' in mol.properties
    assert 'forces' in mol.properties
