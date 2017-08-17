import pytest
import sys


__PYTEST_MARK__ = 'internal'  # mark all tests in this module with this label (see ./conftest.py)


@pytest.mark.tryfirst
def test_lazy_imports():
    import moldesign

    for mod in 'simtk simtk.openmm pyscf pdbfixer scipy':
        assert mod not in sys.modules


def test_package_install_check():
    import moldesign as mdt
    from moldesign.compute import packages

    try:
        from simtk import openmm
    except ImportError:
        assert not packages.openmm.is_installed()
    else:
        assert packages.openmm.is_installed()
        assert packages.openmm.installed_version() == openmm.__version__

