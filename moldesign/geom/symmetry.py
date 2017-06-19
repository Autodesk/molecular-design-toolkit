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

import moldesign as mdt
from moldesign import units as u
from moldesign.external import transformations as trns
from moldesign.interfaces import symmol_interface as smi
from . import __all__ as __toplevel__

get_symmetry = smi.run_symmol
__toplevel__.append('get_symmetry')


class SymmetryElement(object):
    """
    Represents an basic, origin-centered symmetry operation: an identity operation (C1/E),
    inversion center (Ci),
    mirror plane (Cs), rotation axis (C2,C3,...), or
    improper rotation (S4, S6, ...)

    Attributes:
        mol (moldesign.Molecule): molecule this symmetry describes
        idx (int): symmetry index
        symbol (str): Schoenflies symbol
        matrix (np.array): symmetry transformation
        csm (mdt.units.MdtQuantity): mean-squared-distance to this symmetry group (i.e.,
           0 if this symmetry is exact)
        max_diff (mdt.unitsMdtQuantity): maximum distance of any one atom to this symmetry group
    """

    def __init__(self, mol, idx, symbol, matrix, csm, max_diff):
        self.mol = mol
        self.symbol = symbol
        self.matrix = matrix
        self.mol = mol
        self.csm = csm
        self.max_diff = max_diff
        self.idx = idx

    def __str__(self):
        return 'SymmetryElement %s, error=%s' % (self.symbol, self.csm)

    def __repr__(self):
        return '<%s>' % self

    def get_axis(self):
        """
        Returns normal of the plane for Cs or axis of rotation for a Cn
        :param symm_elem: symmetry element (with attributes 'symbol' and 'matrix')
        :return: array of shape (3,)
        """
        mat = np.identity(4)
        mat[:3, :3] = self.matrix
        symbol = self.symbol

        if symbol == 'Cs':
            point, normal = trns.reflection_from_matrix(mat)
            assert np.allclose(point[:3], np.zeros(3))
            return normal

        elif symbol[0] == 'C' and symbol[-1].isdigit():
            angle, normal, point = trns.rotation_from_matrix(mat)
            assert np.allclose(point[:3], np.zeros(3))
            return normal

        elif symbol[0] == 'S' and symbol[-1].isdigit():
            normal = improper_axis_from_matrix(self.matrix)
            return normal

        else:
            raise ValueError('Unrecognized symmetry type %s' % self.symbol)


class MolecularSymmetry(object):
    def __init__(self, mol, symbol, rms,
                 orientation, elems,
                 _job=None):
        self.mol = mol
        self.symbol = symbol
        self.rms = rms
        self.orientation = mdt.utils.if_not_none(orientation, mol.positions)
        self.elems = mdt.utils.if_not_none(elems, [])
        self.groups = mdt.utils.Categorizer(lambda x:x.symbol, self.elems)
        self._job = _job

    @property
    def exact(self):
        """ List[SymmetryElement]: Exact symmetry elements
        """
        return [elem for elem in self.elems if elem.max_diff == 0.0]

    @property
    def approximate(self):
        """ List[SymmetryElement]: Approximate symmetry elements
        """
        return [elem for elem in self.elems if elem.max_diff != 0.0]

    def __str__(self):
        return '%d symmetry element%s' % (len(self.elems), 's' if len(self.elems) != 1 else '')

    def __repr__(self):
        return '%s of molecule %s' % (self, self.mol.name)

    def get_symmetrized_coords(self, elem):
        """
        Symmetrize the molecule based on the symmetry operation
        This will work as long as the symmetry operation brings each atom closest to a symmetry relation.
        """
        import scipy.spatial.distance

        # First, apply the transformation
        oriented_coords = self.orientation
        transformed_coords = self.orientation.T.ldot(elem.matrix).T

        # Next, calculate the correspondence between the untransformed and transformed atoms
        align_to_transform = {}  # map between the original positions and their transformed positions
        transform_to_align = {}  # inverse
        byelement = mdt.utils.Categorizer(lambda x: x.element, self.mol.atoms)
        for elemname, atoms in byelement.items():
            indices = np.array([atom.index for atom in atoms])
            atoms_aligned = oriented_coords[indices].defunits_value()
            atoms_transformed = transformed_coords[indices].defunits_value()
            distances = scipy.spatial.distance.cdist(atoms_aligned, atoms_transformed)
            for a_idx, t_idx in enumerate(distances.argmin(axis=0)):
                align_to_transform[indices[a_idx]] = indices[t_idx]
            for t_idx, a_idx in enumerate(distances.argmin(axis=1)):
                transform_to_align[indices[t_idx]] = indices[a_idx]

        # Make the positions exactly symmetric by averaging them
        pos = np.zeros(transformed_coords.shape) * u.default.length
        for align_atom, transform_atom in align_to_transform.items():
            assert transform_to_align[transform_atom] == align_atom, \
                'Molecule is too far from this symmetry to symmetrize'
            pos[transform_atom] = (oriented_coords[transform_atom] +
                                   transformed_coords[align_atom]) / 2.0
        return pos


def improper_axis_from_matrix(matrix):
    """
    Return rotation angle and axis / mirror plane normal from improper rotation matrix.
    """
    R = np.array(matrix, dtype=np.float64, copy=False)
    R33 = R[:3, :3]
    # direction: unit eigenvector of R33 corresponding to eigenvalue of 1
    w, W = np.linalg.eig(R33.T)
    i = np.where(abs(np.real(w) + 1.0) < 1e-8)[0]
    if not len(i):
        raise ValueError("no unit eigenvector corresponding to eigenvalue -1")
    direction = np.real(W[:, i[-1]]).squeeze()
    return direction
