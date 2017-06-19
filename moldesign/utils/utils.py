from __future__ import print_function, absolute_import, division
from future.builtins import *
from future import standard_library
standard_library.install_aliases()

import future.utils

from functools import reduce
import contextlib
import fractions
import operator
import os
import re
import string
import sys
import tempfile
from html.parser import HTMLParser
from io import StringIO
from uuid import uuid4



def if_not_none(item, default):
    """ Equivalent to `item if item is not None else default` """
    if item is None:
        return default
    else:
        return item


class MLStripper(HTMLParser):
    """ Strips markup language tags from a string.

    FROM http://stackoverflow.com/a/925630/1958900
    """
    def __init__(self):
        if not future.utils.PY2:
            super().__init__()
        self.reset()
        self.fed = []
        self.strict = False
        self.convert_charrefs = True

    def handle_data(self, d):
        self.fed.append(d)

    def get_data(self):
        return ''.join(self.fed)


def html_to_text(html):
    """
    FROM http://stackoverflow.com/a/925630/1958900
    """
    s = MLStripper()
    s.unescape = True  # convert HTML entities to text
    s.feed(html)
    return s.get_data()


def printflush(s, newline=True):
    if newline:
        print(s)
    else:
        print(s, end=' ')
    sys.stdout.flush()


class methodcaller(object):
    """The pickleable implementation of the standard library operator.methodcaller.

    This was copied without modification from:
    https://github.com/python/cpython/blob/065990fa5bd30fb3ca61b90adebc7d8cb3f16b5a/Lib/operator.py

    The c-extension version is not pickleable, so we keep a copy of the pure-python standard library
    code here. See https://bugs.python.org/issue22955

    Original documentation:
    Return a callable object that calls the given method on its operand.
    After f = methodcaller('name'), the call f(r) returns r.name().
    After g = methodcaller('name', 'date', foo=1), the call g(r) returns
    r.name('date', foo=1).
    """
    __slots__ = ('_name', '_args', '_kwargs')

    def __init__(*args, **kwargs):
        if len(args) < 2:
            msg = "methodcaller needs at least one argument, the method name"
            raise TypeError(msg)
        self = args[0]
        self._name = args[1]
        if not isinstance(self._name, future.utils.native_str):
            raise TypeError('method name must be a string')
        self._args = args[2:]
        self._kwargs = kwargs

    def __call__(self, obj):
        return getattr(obj, self._name)(*self._args, **self._kwargs)

    def __repr__(self):
        args = [repr(self._name)]
        args.extend(list(map(repr, self._args)))
        args.extend('%s=%r' % (k, v) for k, v in list(self._kwargs.items()))
        return '%s.%s(%s)' % (self.__class__.__module__,
                              self.__class__.__name__,
                              ', '.join(args))

    def __reduce__(self):
        if not self._kwargs:
            return self.__class__, (self._name,) + self._args
        else:
            from functools import partial
            return partial(self.__class__, self._name, **self._kwargs), self._args


class textnotify(object):
    """ Print a single, immediately flushed line to log the execution of a block.

    Prints 'done' at the end of the line (or 'ERROR' if an uncaught exception)

    Examples:
        >>> import time
        >>> with textnotify('starting to sleep'):
        >>>     time.sleep(3)
        starting to sleep...done
        >>> with textnotify('raising an exception...'):
        >>>     raise ValueError()
        raising an exception...error
        ValueError [...]
    """
    def __init__(self, startmsg):
        if startmsg.strip()[-3:] != '...':
            startmsg = startmsg.strip() + '...'
        self.startmsg = startmsg


    def __enter__(self):
        printflush(self.startmsg, newline=False)

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is None:
            printflush('done')
        else:
            printflush('ERROR')



class BaseTable(object):
    def __init__(self, categories, fileobj=None):
        self.categories = categories
        self.lines = []
        self.fileobj = fileobj

    def add_line(self, obj):
        if hasattr(obj, 'keys'):
            newline = [obj.get(cat, '') for cat in self.categories]
        else:
            assert len(obj) == len(self.categories)
            newline = obj
        self.lines.append(newline)
        self.writeline(newline)

    def writeline(self, newline):
        raise NotImplementedError()

    def getstring(self):
        raise NotImplementedError()


class MarkdownTable(BaseTable):
    def __init__(self, *categories):
        super().__init__(categories)

    def markdown(self, replace=None):
        if replace is None: replace = {}
        outlines = ['| ' + ' | '.join(self.categories) + ' |',
                    '|-' + ''.join('|-' for x in self.categories) + '|']

        for line in self.lines:
            nextline = [str(replace.get(val, val)) for val in line]
            outlines.append('| ' + ' | '.join(nextline) + ' |')
        return '\n'.join(outlines)

    def writeline(self, newline):
        pass

    def getstring(self):
        return self.markdown()


def binomial_coefficient(n, k):
    # credit to http://stackoverflow.com/users/226086/nas-banov
    return int(reduce(operator.mul,
                      (fractions.Fraction(n - i, i + 1) for i in range(k)), 1))


def pairwise_displacements(a):
    """
    :type a: numpy.array
    from http://stackoverflow.com/questions/22390418/pairwise-displacement-vectors-among-set-of-points
    """
    import numpy as np

    n = a.shape[0]
    d = a.shape[1]
    c = binomial_coefficient(n, 2)
    out = np.zeros((c, d))
    l = 0
    r = l + n - 1
    for sl in range(1, n):  # no point1 - point1!
        out[l:r] = a[:n - sl] - a[sl:]
        l = r
        r += n - (sl + 1)
    return out


def is_printable(s):
    import string
    for c in s:
        if c not in string.printable:
            return False
    else:
        return True


class _RedirectStream(object):
    """From python3.4 stdlib
    """

    _stream = None

    def __init__(self, new_target):
        self._new_target = new_target
        # We use a list of old targets to make this CM re-entrant
        self._old_targets = []

    def __enter__(self):
        self._old_targets.append(getattr(sys, self._stream))
        setattr(sys, self._stream, self._new_target)
        return self._new_target

    def __exit__(self, exctype, excinst, exctb):
        setattr(sys, self._stream, self._old_targets.pop())


class redirect_stderr(_RedirectStream):
    """From python3.4 stdlib"""
    _stream = "stderr"


GETFLOAT = re.compile(r'-?\d+(\.\d+)?(e[-+]?\d+)')  # matches numbers, e.g. 1, -2.0, 3.5e50, 0.001e-10


def from_filepath(func, filelike):
    """Run func on a temporary *path* assigned to filelike"""
    if type(filelike) == str:
        return func(filelike)
    else:
        with tempfile.NamedTemporaryFile() as outfile:
            outfile.write(filelike.read().encode())  # hack - prob need to detect bytes
            outfile.flush()
            result = func(outfile.name)
        return result

