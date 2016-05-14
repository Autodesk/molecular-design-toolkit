# Copyright 2016 Autodesk Inc.
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
from copy import deepcopy

from moldesign.units import *


def grid_map(f,v,dims,grids):
    """
    Map function values along a grid
    :param f: function to be evaluated, call signature f(v)
    :param v: vector that sets the static coordinates
    :param dims: ndims-length list of dimensions to vary
    :param grids: ndims-length list of grid values for each dimension
    :return: function value grid
    """
    vmod = deepcopy(v)
    for idx, vals in enumerate(zip(*[g.flat for g in grids])):
        for idim, val in zip(dims, vals): vmod[idim] = val
        if idx == 0:
            firstf = f(vmod)
            gridZ = np.zeros(grids[0].shape) * firstf
            gridZ.flat[0] = firstf
        else:
            gridZ.flat[idx] = f(vmod)
    return gridZ


def function_slice(f,v,dims,ranges):
    """
    Return an arbitrary dimensional slice of function values
    :param f: function to be evaluated, call signature f(v)
    :param v: vector that sets the static coordinates
    :param dims: ndims-length list of dimensions to vary
    :param ranges: ndims-list of values along those dimensions
    :return: gridpoints, function values
    """

    assert len(dims)==len(ranges)
    if len(ranges)>1:
        grids = np.meshgrid(*ranges)
    else:
        grids=list(ranges)

    for igrid,(r,g) in enumerate(zip(ranges,grids)):
        grids[igrid] = units_transfer(r,g)

    gridZ = grid_map(f,v,dims,grids)
    return grids,gridZ

