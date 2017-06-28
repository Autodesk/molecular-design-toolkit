import pytest
import sys


__PYTEST_MARK__ = 'internal'  # mark all tests in this module with this label (see ./conftest.py)


@pytest.mark.tryfirst
def test_lazy_imports():
    import moldesign

    for mod in 'simtk simtk.openmm pyscf pdbfixer scipy':
        assert mod not in sys.modules
