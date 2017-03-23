import pytest

from moldesign import utils


def test_categorizer():
    a = [1,2,3]
    b = ['a', 'b', 'c']
    d = {}

    mycats = utils.Categorizer(len, [a,d])
    mycats.add(b)

    assert set(mycats.keys()) == set([0,3])
    assert len(mycats[0]) == 1
    assert mycats[0][0] is d

    assert len(mycats[3]) == 2
    assert a in mycats[3]
    assert b in mycats[3]


class O(object):
    def __init__(self):
        self.index = None

    def __repr__(self):
        return '<O %s>' % self.index


def test_autoindexlist():

    def test_indices(al):
        for i, o in enumerate(al):
            assert o.index == i

    l = [O(), O(), O()]
    al = utils.AutoIndexList(l)
    test_indices(al)

    oldo = al.pop()
    test_indices(al)

    al.append(O())
    test_indices(al)

    al.append(oldo)
    test_indices(al)

    al.extend([O(), O(), O()])
    test_indices(al)

    al[2] = al.pop(0)
    test_indices(al)

    al[3:5] = [O(), O(), O()]
    test_indices(al)

    al.remove(al[4])
    test_indices(al)

    al.sort(key=id)
    test_indices(al)

    al.insert(0, O())
    test_indices(al)


def test_exclusivelist():
    a = [1,2,3]
    b = ['a', 'b', 'c']
    c = 'hi'
    d = {}

    mylist = utils.ExclusiveList(key=len, iterable=[d, a])

    with pytest.raises(KeyError):
        mylist.append(b)

    with pytest.raises(KeyError):
        mylist.extend([b])

    mylist.remove(a)
    mylist.append(b)

    assert mylist[0] is d
    assert mylist[1] is b
    assert len(mylist) == 2

    with pytest.raises(KeyError):
        mylist.append('abc')

    with pytest.raises(NotImplementedError):
        mylist[0:3] = ['a', 'aa', 'aa']

    with pytest.raises(NotImplementedError):
        mylist.insert(3, 'b')

    newlist = mylist[0:1]
    assert type(mylist) is type(newlist)
    assert mylist._keyfn is newlist._keyfn
    newlist.append(b)
    with pytest.raises(KeyError):
        newlist.append(a)

