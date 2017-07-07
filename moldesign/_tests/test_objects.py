import pickle

import pytest
import numpy as np

from moldesign.utils import Alias
from .object_fixtures import *
from .molecule_fixtures import *
from .object_fixtures import TESTDICT


__PYTEST_MARK__ = 'internal'  # mark all tests in this module with this label (see ./conftest.py)


class ComposedClass(object):
    delegated = Alias('s.lower')

    def __init__(self):
        self.s = 'ABC'


@pytest.mark.parametrize('objkey', pickleable)
def test_pickling(objkey, request):
    obj = request.getfixturevalue(objkey)
    for iprotocol in (0,1,2):
        x = pickle.dumps(obj, protocol=iprotocol)
        y = pickle.loads(x)
        assert type(y) == type(obj)


def test_alias():
    t = ComposedClass()
    assert t.delegated() == 'abc'


def test_dotdict_get(dotdict):
    dd = dotdict
    assert dd == dotdict
    assert dd.get('c', None) == 3
    assert dd.get('ccc', None) is None
    assert dd.a == TESTDICT['a']
    assert dd.d == TESTDICT['d']


def test_dotdict_iterators(dotdict):
    dd = dotdict
    assert len(dd) == len(TESTDICT)
    assert set(dd.keys()) == set(TESTDICT.keys())
    assert set(dd.values()) == set(TESTDICT.values())
    assert dd

    for item in TESTDICT:
        assert item in dd


def test_dotdict_copy(dotdict):
    dd = dotdict.copy()
    assert dd == dotdict
    assert dd._od == dotdict._od

    dd.a = 4
    assert dd.a == dd['a'] == 4
    assert dotdict.a == dotdict['a'] == 1

    del dd.a
    assert 'a' not in dd
    assert 'a' in dotdict


@pytest.mark.screening
def test_dotdict_removals(dotdict):
    dd = dotdict.copy()
    assert dd.pop('_a-a-a', None) is None
    assert dd.pop('d', 4) == 'e'
    assert 'd' not in dd

    dd.d = 5
    assert 'd' in dd

    # item deletion
    del dd.d
    assert 'd' not in dd
    assert len(dd) == len(TESTDICT) - 1
    del dd[3]
    assert 3 not in dd


def test_dotdict_clear(dotdict):
    dd = dotdict.copy()
    dd.clear()
    assert len(dd) == 0
    assert not dd
    assert 'd' not in dd
    assert 'c' not in dd
    assert 'a' not in dd
    with pytest.raises(AttributeError):
        dd.c
    with pytest.raises(KeyError):
        dd['c']


def test_dotdict_introspection(dotdict):
    dd = dotdict.copy()
    assert not hasattr(dd, 'abcd')
    assert hasattr(dd, 'c')
    assert getattr(dd, 'c', None) is 3
    assert getattr(dd, 'ddd', None) is None
    with pytest.raises(AttributeError):
        dd.ddd
    with pytest.raises(KeyError):
        dd['ddd']
    dd['d'] = 'e'


def test_dotdict_consistency(dotdict):
    dd = dotdict.copy()
    dd['k'] = 12345
    assert getattr(dd, 'k') == 12345
    setattr(dd, 'newkey', -42)
    assert dd['newkey'] == -42


def test_dotdict_preserves_ordering(dotdict):
    assert list(dotdict.keys()) == list(TESTDICT.keys())
    assert list(dotdict.values()) == list(TESTDICT.values())
    assert list(dotdict.items()) == list(TESTDICT.items())


@pytest.mark.screening
def test_eigenspace_with_ndarray_identity_permutation():
    from moldesign.mathutils import Eigenspace

    evals = np.arange(3)
    evecs = np.array([[0,1,0],
                      [1,0,0],
                      [0,0,1]],
                     dtype='float')
    espace = Eigenspace(evals, evecs)

    assert (espace.transform([1,0,0]) == [0,1.0, 0.0]).all()
    assert (espace.transform(np.identity(3)) == espace.evecs).all()
