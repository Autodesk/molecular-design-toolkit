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
from math import sqrt
from .. import units as u
from ..utils import exports

@exports
def rmsd_align(src_points, dst_points):
    """ Calculate the 4x4 transformation matrix that aligns a set of source points 
        to a set of destination of points.

        This function implements the method from Kabsch (Acta Cryst. (1978) A34, 827-828) that 
        calculate the optimal rotation matrix that minimizes the root mean squared deviation between 
        two paired sets of points.

    Args:
        src_points(numpy.ndarray): An Nx3 array of the 3D source points.
        dst_points(numpy.ndarray): An Nx3 array of the 3D destination points.

    Returns:
        rmsd_error(Float): The root mean square deviation measuring the alignment error.
        Xform (np.ndarray): The 4x4 array that transforms the source points to the destination points. 
    """

    # Calculate point centers. 
    num_points = src_points.shape[0]
    src_center = src_points.sum(axis=0) / num_points
    dst_center = dst_points.sum(axis=0) / num_points

    # Calculate correlation matrix. 
    C = np.dot(np.transpose(dst_points[:]-dst_center), src_points[:]-src_center)

    # Compute singular value decomposition.
    U, S, V = np.linalg.svd(C)
    d = np.linalg.det(V) * np.linalg.det(U)

    # Make sure the rotation preserves orientation (det = 1). 
    if d < 0.0:
        U[:, -1] = -U[:, -1]

    # Calculate matrix rotating src to dst.
    R = np.dot(U, V)

    # Calculate the 4x4 matrix transforming src to dst.
    Tsrc = np.identity(4, dtype=float)
    Tsrc[0:3,3] = src_center
    Tcenter = np.identity(4, dtype=float)
    Tcenter[0:3,3] = dst_center - src_center
    Rn = np.identity(4, dtype=float)
    Rn[:3,:3] = R
    M1 = np.dot(Rn, np.linalg.inv(Tsrc))
    M2 = np.dot(Tcenter, Tsrc)
    Xform = np.dot(M2, M1)
    print(">>> src center " + str(src_center))
    print(">>> dst center " + str(dst_center))

    # Calculate rmsd error.
    rmsd_error = 0.0
    for i in range(num_points):
        cdiff = R.dot(src_points[i]-src_center) - (dst_points[i]-dst_center)
        rmsd_error += cdiff.dot(cdiff)
    #__for i in range(num_points)
    rmsd_error = sqrt(rmsd_error / num_points)
    return rmsd_error, Xform

@exports
def calculate_rmsd(points1, points2, centered=False):
    """ Calculate the root mean square deviation (RMSD) between two sets of 3D points.

    Args:
        points1(numpy.ndarray): An Nx3 array of the 3D points.
        points2(numpy.ndarray): An Nx3 array of the 3D points.
        centered(Boolean): If true then the points are centered around (0,0,0).

    Returns:
        rmsd(Float): The root mean square deviation. 
    """

    if points1.shape[0] != points2.shape[0]: 
        raise ValueError('The number of points for calculating the RMSD between the first array %d does not equal the number of points in the second array %d' % (num_points1, num_points2))

    # Calculate point centers. 
    num_points = points1.shape[0]
    if centered:
        center1 = np.array([0.0, 0.0, 0.0],dtype=float) 
        center2 = np.array([0.0, 0.0, 0.0],dtype=float) 
    else:
        center1 = points1.sum(axis=0) / num_points
        center2 = points2.sum(axis=0) / num_points

    # Calculate RMSD.
    rmsd = 0.0
    center = center2 - center1 
    for i in range(num_points):
        cdiff = points1[i] - points2[i] + center
        rmsd += cdiff.dot(cdiff)

    return sqrt(rmsd / num_points)

