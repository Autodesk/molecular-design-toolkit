import pytest

import moldesign as mdt

from .molecule_fixtures import small_molecule


registered_types = {}

def typedfixture(*types, **kwargs):
    """This is a decorator that lets us associate fixtures with one or more arbitrary types.
    We'll later use this type to determine what tests to run on the result"""

    def fixture_wrapper(func):
        for t in types:
            registered_types.setdefault(t, []).append(func.__name__)
        return pytest.fixture(**kwargs)(func)

    return fixture_wrapper


@typedfixture('hasmodel')
def openbabel_mmff94(small_molecule):
    small_molecule.set_energy_model(mdt.models.OpenBabelPotential, forcefield='mmff94')
    return small_molecule


@typedfixture('hasmodel')
def openbabel_mmff94s(small_molecule):
    small_molecule.set_energy_model(mdt.models.OpenBabelPotential, forcefield='mmff94s')
    return small_molecule


@typedfixture('hasmodel')
def openbabel_ghemical(small_molecule):
    small_molecule.set_energy_model(mdt.models.OpenBabelPotential, forcefield='ghemical')
    return small_molecule


# This test (along with the uff energy model) is disabled because it does not appear to return a
# gradient that's consistent with the energy surface
@typedfixture('hasmodel')
@pytest.mark.xfail("OpenBabel's UFF implementation appears to return incorrect gradients")
def openbabel_uff(small_molecule):
    small_molecule.set_energy_model(mdt.models.OpenBabelPotential, forcefield='uff')
    return small_molecule


@pytest.mark.parametrize('fixture', registered_types['hasmodel'])
def test_ob_energy_models(request, fixture):
    mol = request.getfixturevalue(fixture)
    assert mol.energy_model is not None
    assert isinstance(mol.calculate_potential_energy(), mdt.units.MdtQuantity)


def test_ob_smiles_read():
    mol = mdt.interfaces.openbabel.from_smiles('CC')
    _assert_its_ethane(mol)


def test_ob_inchi_read():
    mol = mdt.interfaces.openbabel.from_inchi('InChI=1S/C2H6/c1-2/h1-2H3')
    _assert_its_ethane(mol)


def _assert_its_ethane(mol):
    assert mol.num_atoms == 8
    assert len(mol.get_atoms(symbol='H')) == 6
    assert len(mol.get_atoms(symbol='C')) == 2
