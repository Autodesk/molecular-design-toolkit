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
