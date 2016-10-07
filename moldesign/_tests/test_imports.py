import pytest
import sys


def test_lazy_imports():
    import moldesign

    for mod in 'simtk simtk.openmm pyscf pdbfixer scipy':
        assert mod not in sys.modules
