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
import collections
import numpy as np
import scipy.spatial.distance
import ipywidgets as ipy

import moldesign as mdt
from moldesign import units as u
from moldesign.interfaces import symmol_interface as smi
from moldesign.external import transformations as trns

get_symmetry = smi.run_symmol


class SymmetryElement(object):
    """
    Represents an basic, origin-centered symmetry operation: an identity operation (C1/E),
    inversion center (Ci),
    mirror plane (Cs), rotation axis (C2,C3,...), or
    improper rotation (S4, S6, ...)
    """

    def __init__(self, symbol, matrix, **kwargs):
        self.symbol = symbol
        self.matrix = matrix
        for kw, val in kwargs.iteritems():
            setattr(self, kw, val)

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
                 orientation=None,
                 elems=None,
                 **kwargs):
        self.mol = mol
        self.symbol = symbol
        self.rms = rms
        self.orientation = mdt.utils.if_not_none(orientation, mol.atoms.position)
        self.elems = mdt.utils.if_not_none(elems, [])
        for kw, val in kwargs.iteritems():
            setattr(self, kw, val)

    def get_symmetrized_coords(self, elem):
        """
        Symmetrize the molecule based on the symmetry operation
        This will work as long as the symmetry operation brings each atom closest to a symmetry relation.
        """

        # First, apply the transformation
        oriented_coords = self.orientation
        transformed_coords = self.orientation.T.ldot(elem.matrix).T

        # Next, calculate the correspondence between the untransformed and transformed atoms
        align_to_transform = {}  # map between the original positions and their transformed positions
        transform_to_align = {}  # inverse
        byelement = mdt.utils.Categorizer(lambda x: x.element, self.mol.atoms)
        for elemname, atoms in byelement.iteritems():
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
        for align_atom, transform_atom in align_to_transform.iteritems():
            assert transform_to_align[transform_atom] == align_atom, \
                'Molecule is too far from this symmetry to symmetrize'
            pos[transform_atom] = (oriented_coords[transform_atom] +
                                   transformed_coords[align_atom]) / 2.0
        return pos


class Symmetrizer(ipy.Box):
    def __init__(self, mol):
        self._current_shapes = []
        self.mol = mol
        self.tolerance = 0.3 * u.angstrom
        self.original_coords = mol.positions.copy()

        self.showing = ipy.HTML()
        self.viewer = mol.draw3d()
        """:type viewer: moldesign.viewer.GeometryViewer"""

        self.description = ipy.HTML()
        self.symm_selector = ipy.Select()
        self.symm_selector.on_trait_change(self.show_symmetry, name='value')

        self.apply_button = ipy.Button(description='Symmetrize')
        self.apply_button.on_click(self.apply_selected_symmetry)

        self.reset_button = ipy.Button(description='Reset')
        self.reset_button.on_click(self.reset_coords)

        self.apply_all_button = ipy.Button(description='Apply all', padding=10)
        self.apply_all_button.on_click(self.set_highest_symmetry)

        self.tolerance_descrip = ipy.HTML(u'<small>tolerance/\u212B</small>',)
        self.tolerance_chooser = ipy.BoundedFloatText(value=self.tolerance.value_in(u.angstrom),
                                                      min=0.0)

        self.recalculate_button = ipy.Button(description='Recalculate')
        self.recalculate_button.on_click(self.coords_changed)

        self.symm_pane = ipy.VBox([self.description,
                                   self.symm_selector,
                                   ipy.HBox([self.apply_button, self.reset_button]),
                                   self.apply_all_button,
                                   ipy.HBox([self.tolerance_chooser, self.recalculate_button]),
                                   self.tolerance_descrip],
                                  width=325)
        self.coords_changed()

        self.hbox = ipy.HBox([ipy.VBox([self.viewer, self.showing]), self.symm_pane])
        super(Symmetrizer, self).__init__([self.hbox])


    def reset_coords(self, *args):
        self.mol.positions = self.original_coords
        self.viewer.append_frame(positions=self.original_coords)
        self.coords_changed()

    def coords_changed(self, *args):
        self.symm_selector.options = {}
        self.description.value = 'Finding symmetries ...'
        self.tolerance = self.tolerance_chooser.value * u.angstrom
        self.symmetry = get_symmetry(self.mol, tolerance=self.tolerance)
        options = collections.OrderedDict()
        for elem in self.symmetry.elems:
            if elem.max_diff.magnitude != 0.0:
                key = '{0}) {1} (error={2.magnitude:.4f} {2.units})'.format(elem.idx, elem.symbol, elem.max_diff)
            else:
                key = '{0}) {1} (exact)'.format(elem.idx, elem.symbol, elem.max_diff)
            options[key] = elem
        self.symm_selector.options = options
        descrip = 'Highest symmetry group: <b>%s</b><br>' % self.symmetry.symbol
        if self.symmetry.rms.magnitude == 0.0:
            descrip += '(Exact)'
        else:
            descrip += 'RMS Error = {:.03P}'.format(self.symmetry.rms)
        self.description.value = descrip
        self.viewer.append_frame(positions=self.symmetry.orientation)
        self.viewer.center()

    def apply_selected_symmetry(self, *args):
        idx = self.symm_selector.value.idx
        elem = self.symmetry.elems[idx]
        newcoords = self.symmetry.get_symmetrized_coords(elem)
        self.mol.atoms.position = newcoords
        if not np.allclose(newcoords, self.symmetry.orientation, atol=1.0e-10):
            self.viewer.append_frame(positions=newcoords)
            self.coords_changed()

    def show_symmetry(self, *args):
        self.showing.value = ''
        if self._current_shapes:
            for s in self._current_shapes: self.viewer.remove(s, render=False)
            self._current_shapes = []
        if self.symm_selector.value is None:
            self.viewer.render()
            return

        elem = self.symm_selector.value
        symbol = elem.symbol

        self.showing.value = '%s visualization not implemented' % symbol

        if symbol == 'C1':
            self.viewer.render()
            self.showing.value = 'Identity operation'
            return

        elif symbol == 'Ci':
            inversion = self.viewer.draw_sphere(np.zeros(3) * u.angstrom,
                                                color='0x4AB4C4',
                                                radius=0.5 * u.angstrom,
                                                opacity=0.85, render=False)
            self._current_shapes.append(inversion)
            self.showing.value = 'Inversion center'

        elif symbol == 'Cs' or (symbol[0] == 'S' and symbol[1].isdigit()):
            axis = elem.get_axis()
            rad = 2.5 * max(self.symmetry.orientation.max(), 3.0 * u.angstrom)
            plane = self.viewer.draw_circle(np.zeros(3),
                                            axis,
                                            radius=rad,
                                            opacity=0.6,
                                            color='0xAB00FE', render=False)
            self._current_shapes.append(plane)
            self.showing.value = 'Mirror plane (normal = %s)' % axis

        if symbol[0] in 'SC' and symbol[1].isdigit():
            axis = elem.get_axis()
            nrot = int(symbol[1])
            projections = self.symmetry.orientation.dot(axis)
            top = axis * max(3.25 * projections.max(), 3.0*u.angstrom)
            bottom = axis * min(2.5 * projections.min(), -2.5*u.angstrom)
            arrow = self.viewer.draw_arrow(start=bottom, end=top,
                                           color='0x00FE03', render=False, opacity=0.8)
            self._current_shapes.append(arrow)
            if symbol[0] == 'S':
                self.showing.value = '%d-fold improper rotation axis (%s)' % (nrot, axis)
            else:
                self.showing.value = '%d-fold rotation axis (%s)' % (nrot, axis)


        self.viewer.render()


    def set_highest_symmetry(self, *args):
        raise NotImplementedError()



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
