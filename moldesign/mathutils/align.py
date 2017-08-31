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

from math import sqrt
import numpy as np
import moldesign as mdt
from ..utils import exports

@exports
def rmsd_align(src_points, dst_points):
    """ Calculate the 4x4 transformation matrix that aligns a set of source points 
    to a set of destination of points.

    This function implements the method from Kabsch (Acta Cryst. (1978) A34, 827-828) that 
    calculate the optimal rotation matrix that minimizes the root mean squared deviation between 
    two paired sets of points.

    Args:
        src_points (numpy.ndarray): An Nx3 array of the 3D source points.
        dst_points (numpy.ndarray): An Nx3 array of the 3D destination points.

    Returns:
        units.Scalar[length]: The root mean square deviation measuring the alignment error.
        numpy.ndarray: The 4x4 array that transforms the source points to the destination points. 
    """

    if (src_points.shape[0] == 0) or (dst_points.shape[0] == 0):
        n = "Source" if src_points.shape[0] == 0 else "Destination"
        print('WARNING: %s points array for RMSD aligning has zero length.' % n)
        return 0.0, np.identity(4, dtype=float)

    if src_points.shape[0] != dst_points.shape[0]:
        raise ValueError(
            'The number of points for calculating the RMSD between the first array %d \n'
            'does not equal the number of points in the second array %d' % (num_points1, num_points2))

    if mdt.units.get_units(src_points) != mdt.units.get_units(dst_points):
        raise ValueError("The source points units '%s' don't match the destination points units '%s'" % 
            (mdt.units.get_units(src_points), mdt.units.get_units(dst_points)))

    # Numpy matrix/vector multiplication results in a vector with dimensionless units 
    # even if the vector has units. To retain the vector's units multiple the result by 
    # the points units.
    if isinstance(src_points, mdt.units.MdtQuantity):
        units = mdt.units.get_units(src_points)
    else:
        units = 1.0

    # Calculate point centers. 
    num_points = src_points.shape[0]
    src_center = src_points.sum(axis=0) / num_points
    dst_center = dst_points.sum(axis=0) / num_points

    # Calculate correlation matrix. 
    corr_mat = np.dot(np.transpose(dst_points[:]-dst_center), src_points[:]-src_center)

    # Compute singular value decomposition.
    u, s, v = np.linalg.svd(corr_mat)
    det = np.linalg.det(v) * np.linalg.det(u)

    # Make sure the rotation preserves orientation (det = 1). 
    if det < 0.0:
        u[:, -1] = -u[:, -1]

    # Calculate matrix rotating src to dst.
    rot_mat = np.dot(u, v)

    # Calculate the 4x4 matrix transforming src to dst.
    tsrc = np.identity(4, dtype=float)
    tsrc[0:3,3] = src_center
    tcenter = np.identity(4, dtype=float)
    tcenter[0:3,3] = dst_center - src_center
    rn = np.identity(4, dtype=float)
    rn[:3,:3] = rot_mat
    m1 = np.dot(rn, np.linalg.inv(tsrc))
    m2 = np.dot(tcenter, tsrc)
    xform = np.dot(m2, m1)

    # Calculate rmsd error.
    rmsd_error = 0.0
    for i in range(num_points):
        cdiff = rot_mat.dot(src_points[i]-src_center)*units - (dst_points[i]-dst_center)
        rmsd_error += cdiff.dot(cdiff)
    #__for i in range(num_points)
    rmsd_error = np.sqrt(rmsd_error / num_points)
    return rmsd_error, xform

@exports
def calculate_rmsd(points1, points2, centered=False):
    """ Calculate the root mean square deviation (RMSD) between two sets of 3D points.

    Args:
        points1 (numpy.ndarray): An Nx3 array of the 3D points.
        points2 (numpy.ndarray): An Nx3 array of the 3D points.
        centered (bool): If true then the points are centered around (0,0,0).

    Returns:
        units.Scalar[length]: The root mean square deviation. 
    """
    if (points1.shape[0] == 0) or (points2.shape[0] == 0):
        n = 1 if points1.shape[0] == 0 else 2
        print('WARNING: Points%d array for RMSD calculation has zero length.' % n)
        return 0.0 

    if points1.shape[0] != points2.shape[0]: 
        raise ValueError(
            'The number of points for calculating the RMSD between the first array %d \n'
            'does not equal the number of points in the second array %d' % (num_points1, num_points2))

    if mdt.units.get_units(points1) != mdt.units.get_units(points2):
        raise ValueError("Points1 units '%s' don't match the points2 units '%s'" %
            (mdt.units.get_units(points1), mdt.units.get_units(points2)))

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

    return np.sqrt(rmsd / num_points)

