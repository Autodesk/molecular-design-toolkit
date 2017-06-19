import pytest
import sys


@pytest.mark.tryfirst
def test_lazy_imports():
    import moldesign

    for mod in 'simtk simtk.openmm pyscf pdbfixer scipy':
        assert mod not in sys.modules
