import pytest

import moldesign as mdt


registered_types = {}


def typedfixture(*types, **kwargs):
    """This is a decorator that lets us associate fixtures with one or more arbitrary types.
    We'll later use this type to determine what tests to run on the result"""

    def fixture_wrapper(func):
        for t in types:
            registered_types.setdefault(t, []).append(func.__name__)
        return pytest.fixture(**kwargs)(func)

    return fixture_wrapper


@pytest.fixture
def ammonium_nocharge():
    return mdt.from_smiles('[NH4]')


@pytest.fixture
def ammonium_charged():
    return mdt.from_smiles('[NH4+]')


@pytest.mark.parametrize('objkey',
                         [ammonium_nocharge, ammonium_charged])
def test_ammonium_formal_charge(objkey, request):
    mol = request.getfuncargvalue(objkey)
    mdt.assign_formal_charges(mol)
    assert mol.charge == 1

    for atom in mol:
        if atom.atnum == 6:
            assert atom.formal_charge == 1
        else:
            assert atom.atnum == 1
            assert atom.formal_charge == 0
