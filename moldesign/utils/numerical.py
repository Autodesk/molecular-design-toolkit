from __future__ import print_function, absolute_import, division
from future.builtins import *
from future import standard_library
standard_library.install_aliases()

# Copyright 2017 Autodesk Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import numpy as np

from . import Alias


class ResizableArray(object):
    """ Behaves like a numpy array, but with fast extends and appends for the first index.

    Similar to python lists, the underlying array is reallocated in large chunks whenever the
    number of elements grows beyond the current memory allocation.
    """
    #@args_from(np.array)
    def __init__(self, *args, **kwargs):
        self._array = np.array(*args, **kwargs)
        self._len = len(self._array)
        self._subarray = self._array[:self._len]
        self._size = self._len

    def __getattr__(self, item):
        if item == '_subarray':
            return self.__getattribute__('_subarray')
        return getattr(self._subarray, item)

    def __repr__(self):
        return '%s(%s)' % (self.__class__.__name__, repr(self._subarray))

    def append(self, item):
        self.extend([item])

    def extend(self, its):
        try:
            ll = len(its)
        except TypeError:
            its = list(its)
            ll = len(its)

        newlen = self._len + ll
        if newlen > self._size:
            self._resize(newlen)

        for item in its:
            self._array[self._len] = item
            self._len += 1

        self._subarray = self._array[:self._len]

    def _resize(self, minsize):
        newsize = round_up_to_power_of_two(minsize)
        if newsize <= self._size:
            return

        newshape = (newsize,) + self._array.shape[1:]
        newarray = np.empty(newshape, dtype=self._array.dtype)
        newarray[:self._len] = self._array

        self._array = newarray
        self._size = newsize

# delegate magic methods as well - this is a list of all math-related array methods
_ARRAYMAGIC = ('__abs__', '__add__', '__and__', '__array__', '__contains__', '__copy__',
               '__deepcopy__', '__delitem__', '__delslice__', '__div__', '__divmod__', '__eq__',
               '__float__', '__floordiv__', '__ge__', '__getitem__', '__getslice__', '__gt__',
               '__hex__', '__iadd__', '__iand__', '__idiv__', '__ifloordiv__', '__ilshift__',
               '__imod__', '__imul__', '__index__', '__int__', '__invert__', '__ior__', '__ipow__',
               '__irshift__', '__isub__', '__itruediv__', '__ixor__', '__le__', '__len__',
               '__long__', '__lshift__', '__lt__', '__mod__', '__mul__', '__ne__', '__neg__',
               '__nonzero__', '__oct__', '__or__', '__pos__', '__pow__', '__radd__', '__rand__',
               '__rdiv__', '__rdivmod__', '__rfloordiv__', '__rlshift__', '__rmod__', '__rmul__',
               '__ror__', '__rpow__', '__rrshift__', '__rshift__', '__rsub__', '__rtruediv__',
               '__rxor__', '__setitem__', '__setslice__', '__str__', '__sub__',
               '__truediv__', '__xor__')

for _methname in _ARRAYMAGIC:
    setattr(ResizableArray, _methname, Alias('_subarray.%s' % _methname))



def round_up_to_power_of_two(n):
    """
    From http://stackoverflow.com/a/14267825/1958900
    """
    if n < 0:
        raise TypeError("Nonzero positive integers only")
    elif n == 0:
        return 1
    else:
        return 1 << (n-1).bit_length()
