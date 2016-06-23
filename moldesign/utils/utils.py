"""
My standard utilities. Intended to be included in all projects
Obviously everything included here needs to be in the standard library (or numpy)
"""
import contextlib
import fractions
import operator
import os
import re
import shutil
import string
import sys
import tempfile
import threading
from cStringIO import StringIO
from uuid import uuid4

import webcolors


def make_none(): return None


@contextlib.contextmanager
def recursionlimit_atleast(n=1000):
    """Context manager for temporarily raising the context manager's
    the interpreter's maximum call stack size (misleading called the ``recursion limit``)

    Notes:
        This will explicitly reset the the recursion limit when we exit the context;
            any intermediate recursion limit changes will be lost
        This will not lower the limit ``n`` is less than the current recursion limit.
    """
    current_limit = sys.getrecursionlimit()
    if n >= current_limit:
        sys.setrecursionlimit(n)
    yield
    sys.setrecursionlimit(current_limit)


def if_not_none(item, default):
    """ Equivalent to `item if item is not None else default` """
    if item is None:
        return default
    else:
        return item


class PipedFile(object):
    """
    Allows us to pass data by filesystem path without ever writing it to disk
    To prevent deadlock, we spawn a thread to write to the pipe
    Call it as a context manager:
    >>> with PipedFile('file contents',filename='contents.txt') as pipepath:
    >>>     print open(pipepath,'r').read()
    """
    def __init__(self, fileobj, filename='pipe'):
        if type(fileobj) in (unicode,str):
            self.fileobj = StringIO(fileobj)
        else:
            self.fileobj = fileobj
        self.tempdir = None
        assert '/' not in filename,"Filename must not include directory"
        self.filename = filename

    def __enter__(self):
        self.tempdir = tempfile.mkdtemp()
        self.pipe_path = os.path.join(self.tempdir, self.filename)
        os.mkfifo(self.pipe_path)
        self.pipe_thread = threading.Thread(target=self._write_to_pipe)
        self.pipe_thread.start()
        return self.pipe_path

    def _write_to_pipe(self):
        with open(self.pipe_path,'w') as pipe:
            pipe.write(self.fileobj.read())

    def __exit__(self, type, value, traceback):
        if self.tempdir is not None:
            shutil.rmtree(self.tempdir)


def remove_directories(list_of_paths):
    """
    Removes non-leafs from a list of directory paths
    """
    found_dirs = set('/')
    for path in list_of_paths:
        dirs = path.strip().split('/')
        for i in xrange(2, len(dirs)):
            found_dirs.add('/'.join(dirs[:i]))

    paths = [path for path in list_of_paths if
             (path.strip() not in found_dirs) and path.strip()[-1] != '/']
    return paths


def make_local_temp_dir():
    tempdir = '/tmp/%s' % uuid4()
    os.mkdir(tempdir)
    return tempdir


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


class PrintTable(BaseTable):
    def __init__(self, formatstr, fileobj=sys.stdout):
        self.format = formatstr
        categories = []
        self._wrote_header = False
        for field in string.Formatter().parse(formatstr):
            key = field.split('.')[0]
            categories.append(key)
        super(PrintTable, self).__init__(categories, fileobj=fileobj)

    def writeline(self, line):
        if not self._wrote_header:
            print >> self._fileobj, self.format.format(self.categories)
            self._wrote_header = True

        if self.fileobj is None: return
        print >> self.fileobj, self.formatstr.format(**line)

    def getstring(self):
        s = StringIO()
        for line in self.lines:
            print >> s, self.format.format(line)
        return s.getvalue()


class MarkdownTable(BaseTable):
    def __init__(self, *categories):
        super(MarkdownTable, self).__init__(categories)

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


class redirect_stdout(_RedirectStream):
    """From python3.4 stdlib"""
    _stream = "stdout"


class redirect_stderr(_RedirectStream):
    """From python3.4 stdlib"""
    _stream = "stderr"


GETFLOAT = re.compile(r'-?\d+(\.\d+)?(e[-+]?\d+)')  # matches numbers, e.g. 1, -2.0, 3.5e50, 0.001e-10


def is_color(s):
    """ Do our best to determine if "s" is a color spec that can be converted to hex
    :param s:
    :return:
    """
    def in_range(i): return 0 <= i <= int('0xFFFFFF', 0)

    try:
        if type(s) == int:
            return in_range(s)
        elif type(s) not in (str, unicode):
            return False
        elif s in webcolors.css3_names_to_hex:
            return True
        elif s[0] == '#':
            return in_range(int('0x' + s[1:], 0))
        elif s[0:2] == '0x':
            return in_range(int(s, 0))
        elif len(s) == 6:
            return in_range(int('0x' + s, 0))
    except ValueError:
        return False