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
from moldesign import units as u
from moldesign.mathutils import perpendicular

from . import GeometryViewer


class BondClicker(GeometryViewer):
    """
    Allow the user to highlight bonds - this is a hack around 3dmol.js to allow
    clickable bonds.
    """
    BONDCOLOR = '#C8C8C8'
    ATOMRADIUS = 0.35
    SINGLERADIUS = 0.18
    DOUBLERADIUS = 0.13
    TRIPLERADIUS = 0.12
    DOUBLEOFFSET = 0.16
    TRIPLEOFFSET = 0.14

    def __init__(self, mol, **kwargs):
        self._bonds = {}
        self._bond_shapes = {}
        self._bond_colors = {}
        super(BondClicker, self).__init__(mol=mol, render=False, **kwargs)
        self.atom_callbacks = []
        self.bond_callbacks = []
        self.click_callbacks = []
        self.vdw(radius=self.ATOMRADIUS, render=False)
        self.draw_all_bonds(render=True)

    def set_positions(self, *args, **kwargs):
        render = kwargs.get('render',True)
        kwargs['render'] = False
        super(BondClicker, self).set_positions(*args, **kwargs)
        self.draw_all_bonds(render=render)

    def draw_all_bonds(self, render=True, batch=True):
        # TODO: this should be written in javascript, for speed and consistency
        for bond in self.mol.bonds:
            self.draw_bond(bond, batch=batch, render=False)
        if render: self.render()

    def set_bond_color(self, color, bond, render=True):
        self._bond_colors[bond] = color
        self.draw_bond(bond, render=render)

    def unset_bond_color(self, bond, render=True):
        self._bond_colors.pop(bond, None)
        self.draw_bond(bond, render=render)

    def draw_bond(self, bond, render=True, batch=False, **shape_args):
        atom = bond.a1
        nbr = bond.a2
        order = bond.order
        assert atom.index != nbr.index

        if bond in self._bond_shapes:  # i.e., we need to remove and redraw this bond
            for shape in self._bond_shapes[bond]:
                self.remove(shape, render=False, batch=batch)

        assert 'clickable' not in shape_args
        color = self._bond_colors.get(bond, self.BONDCOLOR)
        kwargs = dict(render=False,
                      clickable=True,
                      color=color,
                      batch=batch)
        kwargs.update(shape_args)

        vec = atom.position - nbr.position
        if order == 2:
            kwargs['radius'] = self.DOUBLERADIUS
            offset = self.DOUBLEOFFSET * perpendicular(vec) * u.angstrom
            t1 = self.draw_tube(atom.position + offset,
                                nbr.position + offset,
                                **kwargs)
            t2 = self.draw_tube(atom.position - offset,
                                nbr.position - offset,
                                **kwargs)
            shapes = [t1, t2]
        elif order == 3:
            kwargs['radius'] = self.TRIPLERADIUS
            offset = self.TRIPLEOFFSET * perpendicular(vec) * u.angstrom
            t1 = self.draw_tube(atom.position + offset,
                                nbr.position + offset,
                                **kwargs)
            t2 = self.draw_tube(atom.position - offset,
                                nbr.position - offset,
                                **kwargs)
            t3 = self.draw_tube(atom.position, nbr.position,
                                **kwargs)
            shapes = [t1, t2, t3]
        else:  # assume single bond
            kwargs['radius'] = self.SINGLERADIUS
            tube = self.draw_tube(atom.position, nbr.position,
                                  **kwargs)
            shapes = [tube]

        self._bond_shapes[bond] = shapes
        for x in shapes:
            self._bonds[x] = bond
        if render: self.render()

    def handle_click(self, trait_name, old, new):
        if 'pyid' in new:
            pyid = new['pyid']
            obj = self._bonds[pyid]
            for func in self.bond_callbacks:
                func(obj)
        else:
            obj = self.mol.atoms[new['index']]
            for func in self.atom_callbacks:
                func(obj)

        for func in self.click_callbacks:
            func(obj)