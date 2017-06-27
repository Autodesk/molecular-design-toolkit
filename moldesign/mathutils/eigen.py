from __future__ import print_function, absolute_import, division
from future import standard_library
standard_library.install_aliases()
from future.builtins import *

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

from .. import units as u
from ..utils import exports
from ..mathutils import normalized


@exports
class Eigenspace(object):
    """ Holds sets of eigenvactors and eigenvalues, offers helpful methods for organizing them.

    Args:
        evals (List[Scalar]): list of eigenvalues
        evecs (List[Vector]): list of eigenvectors (in the same order as eigenvalues). These will
           automatically be normalized

    Note:
        The eigenvectors from many eigenvector solvers - notably including scipy's - will need to be
        transposed to fit the form of ``evecs`` here!

        It's generally assumed that Eigenspace objects will be subclassed to wrap various bits of
        eignproblem-related functionality.
    """
    def __init__(self, evals, evecs):
        if len(evals) != len(evecs):
            raise ValueError('Eigenvalues and eigenvectors have different lengths!')

        self.evals = u.array(evals)
        self.evecs = u.array([normalized(evec, zero_as_zero=True)
                              for evec in evecs])

    def __str__(self):
        return "%s of dimension %s" % (self.__class__.__name__, len(self.evals))

    def sort(self, largest_first=True, key=abs):
        """ Sort the eigenvectors and values in place*

        By default, this sorts from largest to smallest by the _absolute magnitude_ of the
        eigenvalues.

        Note:
            *this sort is only "in place" in the sense that it mutates the data in this instance;
              note that it still uses auxiliary memroy

        Args:
            largest_first (bool): sort from largest to smallest (equivalent to ``reverse=True`` in a
                standard python sort, except that it is true by default here)
            key (callable): function of the form ``f(eigenval)`` OR ``f(eval, evec)``.
                  By default, sorts by the ``abs`` of the eigenvalues
        """
        try:  # construct the sorting function
            result = key(self.evals[0], self.evecs[0])
        except TypeError:
            def keyfn(t):
                return key(t[0])
        else:
            def keyfn(t):
                return key(t[0], t[1])

        evals, evecs = zip(*sorted(zip(self.evals.copy(), self.evecs.copy()),
                                   key=keyfn, reverse=largest_first))
        self.evals[:] = evals
        self.evecs[:] = evecs

    def transform(self, coords):
        """ Transform a list of coordinates (or just a single one) into this eigenbasis

        Args:
            coords (List[Vector] or Vector): coordinate(s) to transform

        Returns:
            List[Vector] or vector: transformed coordinate(s)
        """
        c = u.array(coords)
        dims = len(c.shape)
        if dims == 1:
            return u.dot(self.evecs, c)
        elif dims == 2:
            return u.dot(self.evecs, c.T).T
        else:
            raise ValueError('Transform accepts only 1- or 2-dimensional arrays')


