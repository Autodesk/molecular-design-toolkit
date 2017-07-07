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
import itertools

import numpy as np

from .. import utils
from .. import units as u


@utils.exports
class VolumetricGrid(object):
    """ Creates a 3D, rectangular grid of points

    Args:
        xrange (Tuple[len=2]): (min,max) in x direction
        yrange (Tuple[len=2]): (min,max) in y direction
        zrange (Tuple[len=2]): (min,max) in z direction
        xpoints (int): number of grid lines in x direction (default: 25)
        ypoints (int): number of grid lines in y direction (default: xpoints)
        zpoints (int): number of grid lines in z direction (default: xpoints)
    """
    dx, dy, dz = (utils.IndexView('deltas', i) for i in range(3))
    xr, yr, zr = (utils.IndexView('ranges', i) for i in range(3))
    xpoints, ypoints, zpoints = (utils.IndexView('points', i) for i in range(3))
    xspace, yspace, zspace = (utils.IndexView('spaces', i) for i in range(3))

    def __init__(self, xrange, yrange, zrange,
                 xpoints=25, ypoints=None, zpoints=None):

        self.points = np.array([xpoints,
                                ypoints if ypoints is not None else xpoints,
                                zpoints if zpoints is not None else xpoints],
                               dtype='int')

        self.ranges = u.array([xrange, yrange, zrange])
        self.deltas = (self.ranges[:,1] - self.ranges[:,0]) / (self.points - 1.0)
        self.spaces = [u.linspace(*r, num=p) for r,p in zip(self.ranges, self.points)]

    @property
    def origin(self):
        """ Vector[len=3]: the origin of the grid (the lowermost corner in each dimension)
        """
        origin = [r[0] for r in (self.ranges)]
        try:
            return u.array(origin)
        except u.DimensionalityError:
            return origin

    @property
    def npoints(self):
        """ int: total number of grid points in this grid
        """
        return np.product(self.points)

    def iter_points(self):
        """ Iterate through every point on the grid.

        Always iterates in the same order. The x-index is the most slowly varying,
        the z-index is the fastest.

        Yields:
            Vector[len=3]: x,y,z coordinate of each point on the grid
        """
        for i,j,k in itertools.product(*map(range, self.points)):
                yield self.origin + self.deltas * [i,j,k]

    def allpoints(self, dtype='float'):
        """ Return an array of all coordinates on the grid.

        This obviously takes a lot of memory, but is useful for evaluating vectorized functions
        on this grid. Points are returned in the same order as ``iter_points``.

        Yields:
            Matrix[shape=(self.npoints**3,3)]: x,y,z coordinate of each point on the grid
        """
        ap = np.empty((self.npoints, 3), dtype=dtype)
        if hasattr(self.xr, 'units'):
            ap = ap * self.xr.units

        for ip, point in enumerate(self.iter_points()):
            ap[ip] = point

        return ap


@utils.exports
def padded_grid(positions, padding, npoints=25):
    """ Creates a 3D, rectangular grid of points surrounding a set of positions
    The points in the grid are evenly spaced in each dimension, but spacing may differ between
    in different dimensions.

    Args:
        positions (Matrix[shape=(*,3)]): positions to create the grid around
        padding (Scalar): how far to extend the grid past the positions in each dimension
        npoints (int): number of points in each direction (total number of points is npoints**3)

    Returns:
        VolumetricGrid: grid object
    """
    mins = positions.min(axis=0)-padding
    maxes = positions.max(axis=0)+padding

    xr = (mins[0], maxes[0])
    yr = (mins[1], maxes[1])
    zr = (mins[2], maxes[2])

    return VolumetricGrid(xr, yr, zr, npoints)
