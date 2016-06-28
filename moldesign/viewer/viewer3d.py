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
import sys

import IPython.display as dsp
import numpy as np

import moldesign as mdt
from moldesign import units as u
from moldesign import utils
from moldesign.helpers import colormap, VolumetricGrid
from moldesign.utils import is_color
from nbmolviz.drivers3d import MolViz_3DMol
from . import toplevel


# Right now hard-coded to use 3DMol driver - will need to be configurable when there's more than one
# In theory, the only tie to the specific driver should be the class that we inherit from
@toplevel
class GeometryViewer(MolViz_3DMol):
    """
    Viewer for static and multiple-frame geometries

    :ivar mol: Buckyball molecule
    :type mol: moldesign.Molecule.molecule
    """
    DISTANCE_UNITS = u.angstrom
    HIGHLIGHT_COLOR = '#1FF3FE'
    DEFAULT_COLOR_MAP = colormap
    DEFAULT_WIDTH = 625
    DEFAULT_HEIGHT = 400

    def __reduce__(self):
        """prevent these from being pickled for now"""
        return utils.make_none, tuple()

    def __init__(self, mol=None, style=None, display=False, render=True, **kwargs):
        """
        TODO: make clickable atoms work with trajectories
        TODO: coloring methods

        Args:
            mol (moldesign.AtomContainer): Atoms to visualize
        """
        kwargs.setdefault('height', self.DEFAULT_HEIGHT)
        kwargs.setdefault('width', self.DEFAULT_WIDTH)
        super(GeometryViewer, self).__init__(**kwargs)
        self.set_click_callback(callback=self.handle_click)
        self.selection_group = None
        self.selection_id = None
        self.atom_highlights = []
        self.frame_change_callback = None
        self.wfn = None
        self._cached_orbitals = set()
        self._callbacks = set()
        self._axis_objects = None
        self._frame_positions = []
        if mol:
            self.add_molecule(mol, render=False)
            self._frame_positions.append(self.get_positions())
            if style is None:
                self.autostyle(render=render)
            else:
                self.set_style(style, render=render)
        if display: dsp.display(self)

    def autostyle(self, render=True):
        if self.mol.mass <= 500.0 * u.dalton:
            self.stick()
        else:
            cartoon_atoms = []
            line_atoms = []
            stick_atoms = []
            biochains = set()
            for residue in self.mol.residues:
                if residue.type in ('protein', 'dna', 'rna'):
                    biochains.add(residue.chain)

                if residue.type == 'protein':
                    cartoon_atoms.extend(residue.atoms)
                elif residue.type in ('water', 'solvent'):
                    line_atoms.extend(residue.atoms)
                elif residue.type in ('dna', 'rna') and self.mol.numatoms > 1000:
                    cartoon_atoms.extend(residue.atoms)
                else:  # includes DNA, RNA if molecule is small enough
                    stick_atoms.extend(residue.atoms)

            if cartoon_atoms:
                self.cartoon(atoms=cartoon_atoms, render=False)
                if len(biochains) > 1:
                    self.color_by('chain', atoms=cartoon_atoms, render=False)
                else:
                    self.color_by('residue.resname', atoms=cartoon_atoms, render=False)
            if line_atoms:
                self.line(atoms=line_atoms, render=False)
            if stick_atoms:
                self.stick(atoms=stick_atoms, render=False)

        # Deal with unbonded atoms (they only show up in VDW rep)
        if self.mol.numatoms > 1000:
            print 'WARN: large structure; waters not shown by default.'
            lone = [atom for atom in self.mol.atoms if
                    atom.num_bonds == 0 and atom.residue.type != 'water']
            self.hide(atoms=[atom for atom in self.mol.atoms if atom.num_bonds == 0 and
                             atom.residue.type == 'water'])
        else:
            self.show_unbonded()

        if render: self.render()

    def show_unbonded(self, radius=0.5):
        lone = [atom for atom in self.mol.atoms if atom.num_bonds == 0]
        if lone: self.vdw(atoms=lone, render=False, radius=radius)

    @staticmethod
    def _atoms_to_json(atomlist):
        if hasattr(atomlist, 'iteratoms'):
            idxes = [a.index for a in atomlist.iteratoms()]
        else:
            # TODO: verify that these are actually atoms and not something else with a .index
            idxes = [a.index for a in atomlist]
        atomsel = {'index': idxes}
        return atomsel

    def _set_bonds(self):
        """
        Sends bond structure to the drawing software
        This is necessary for formats like PDB, where the bond structure is set implicitly
        """
        bonds = []
        for atom in self.mol.atoms:
            nbrs = atom.bond_graph.keys()
            bonds.append({'index': atom.index,
                          'nbr': [n.index for n in nbrs],
                          'order': [atom.bond_graph[n] for n in nbrs]})
        self.viewer('setBonds', [bonds])

    def color_by(self, atom_callback, atoms=None, mplmap='auto', render=True):
        """
        Color atoms according to either:
          * an atomic attribute (e.g., 'chain', 'residue', 'mass')
          * a callback function that accepts an atom and returns a color or a category
        :param atom_callback: callable f(atom) returns color or category
        :returns: dict of categories
        """
        atoms = utils.if_not_none(atoms, self.mol.atoms)
        if isinstance(atom_callback, str):  # shortcut to use strings to access atom attributes, i.e. "ff.partial_charge"
            attrs = atom_callback.split('.')

            def atom_callback(atom):
                obj = atom
                for attr in attrs:
                    obj = getattr(obj, attr)
                return obj

        if callable(atom_callback):
            colors = utils.Categorizer(atom_callback, atoms)

        name_is_color = map(is_color, colors.keys())
        if len(colors) <= 1:
            colors = {'gray': atoms}

        elif not all(name_is_color):
            assert not any(name_is_color),\
                "callback function returned a mix of colors and categories"
            categories = colors
            cats = categories.keys()
            # If there are >256 categories, this is a many-to-one mapping
            colornames = colormap(cats, mplmap=mplmap)
            colors = {c: [] for c in colornames}
            for cat, color in zip(cats, colornames):
                colors[color].extend(categories[cat])

        self.set_colors(colors, render=render)

    def get_input_file(self):
        if len(self.mol.atoms) <= 250:
            fmt = 'sdf'
        else:
            fmt = 'pdb'
        if not hasattr(self.mol, 'write'):
            writemol = mdt.Molecule(self.mol)
        else:
            writemol = self.mol
        instring = writemol.write(format=fmt)
        return instring, fmt

    def get_positions(self):
        positions = [atom.position.value_in(self.DISTANCE_UNITS).tolist()
                     for atom in self.mol.atoms]
        return positions

    def calc_orb_grid(self, orbname, npts=30, padding=2.5*u.angstrom):
        """ Calculate orbitals on a grid

        Args:
            orbname: Either orbital index (for canonical orbitals) or a tuple ( [orbital type],
               [orbital index] ) where [orbital type] is a keyword (e.g., canonical, natural, nbo,
               ao, etc)
            npts (int): number of points in each dimension of the grid
            padding (u.Scalar[length]): extent of the grid beyond the furthest atom in each
                dimension

        Returns:
            moldesign.viewer.VolumetricGrid: orbital values on a grid
        """
        # TODO: why is orbname sometimes an np.int64 instead of an int?
        # NEWFEATURE: limit bounding box based on the non-zero atomic centers. Useful for localized orbitals,
        #             which otherwise require high resolution
        try:
            orbtype, orbidx = orbname
        except TypeError:
            orbtype = 'canonical'
            orbidx = orbname

        # self.wfn should already be set to the wfn for the current frame
        orbital = self.wfn.orbitals[orbtype][orbidx]
        positions = self._frame_positions[self.current_frame]*u.angstrom
        grid = VolumetricGrid(positions,
                              padding=padding, npoints=npts)
        coords = grid.xyzlist().reshape(3, grid.npoints ** 3).T
        values = orbital(coords)
        grid.fxyz = values.reshape(npts, npts, npts)
        return grid

    def get_orbnames(self):
        raise NotImplementedError()

    def draw_atom_vectors(self, vecs, rescale_to=1.75,
                          scale_factor=None, opacity=0.85,
                          radius=0.11, render=True, **kwargs):
        """
        For displaying atom-centered vector data (e.g., momenta, forces)
        :param rescale_to: rescale to this length (in angstroms) (not used if scale_factor is passed)
        :param scale_factor: Scaling factor for arrows: dimensions of [vecs dimensions] / [length]
        :param render: render immediately
        :param kwargs: keyword arguments for self.draw_arrow
        """
        kwargs['radius'] = radius
        kwargs['opacity'] = opacity
        if vecs.shape == (self.mol.ndims,):
            vecs = vecs.reshape(self.mol.num_atoms, 3)
        assert vecs.shape == (self.mol.num_atoms, 3), '`vecs` must be a num_atoms X 3 matrix or' \
                                                      ' 3*num_atoms vector'
        assert not np.allclose(vecs, 0.0), "Input vectors are 0."

        # strip units and scale the vectors appropriately
        if scale_factor is not None:  # scale all arrows by this quantity
            if (u.get_units(vecs)/scale_factor).dimensionless:  # allow implicit scale factor
                scale_factor = scale_factor / self.DISTANCE_UNITS

            vecarray = vecs / scale_factor
            try: arrowvecs = vecarray.value_in(self.DISTANCE_UNITS)
            except AttributeError: arrowvecs = vecarray

        else:  # rescale the maximum length arrow length to rescale_to
            try:
                vecarray = vecs.magnitude
                unit = vecs.getunits()
            except AttributeError:
                vecarray = vecs
                unit = ''
            lengths = np.sqrt((vecarray * vecarray).sum(axis=1))
            scale = (lengths.max() / rescale_to)  # units of [vec units] / angstrom
            if hasattr(scale,'defunits'): scale = scale.defunits()
            arrowvecs = vecarray / scale
            print 'Arrow scale: {q:.3f} {unit} per {native}'.format(q=scale, unit=unit,
                                                                    native=self.DISTANCE_UNITS)
        shapes = []
        for atom, vecarray in zip(self.mol.atoms, arrowvecs):
            shapes.append(self.draw_arrow(atom.position, vector=vecarray, render=False, **kwargs))
        if render: self.render()
        return shapes

    def draw_axis(self, on=True, render=True):
        label_kwargs = dict(color='white', opacity=0.4, render=False, fontsize=14)
        if on and self._axis_objects is None:
            xarrow = self.draw_arrow([0, 0, 0], [1, 0, 0], color='red', render=False)
            xlabel = self.draw_label([1.0, 0.0, 0.0], text='x', **label_kwargs)
            yarrow = self.draw_arrow([0, 0, 0], [0, 1, 0], color='green', render=False)
            ylabel = self.draw_label([-0.2, 1, -0.2], text='y', **label_kwargs)
            zarrow = self.draw_arrow([0, 0, 0], [0, 0, 1], color='blue', render=False)
            zlabel = self.draw_label([0, 0, 1], text='z', **label_kwargs)
            self._axis_objects = [xarrow, yarrow, zarrow,
                                  xlabel, ylabel, zlabel]

        elif not on and self._axis_objects is not None:
            for arrow in self._axis_objects:
                self.remove(arrow, render=False)
            self._axis_objects = None

        if render: self.render()

    def draw_forces(self, **kwargs):
        return self.draw_atom_vectors(self.mol.forces, **kwargs)

    def draw_momenta(self, **kwargs):
        return self.draw_atom_vectors(self.mol.momenta, **kwargs)

    def highlight_atoms(self, atoms, render=True):
        # TODO: Need to handle style changes
        if self.atom_highlights:
            self.unset_color(self.atom_highlights)
        self.atom_highlights = []
        self.set_color(self.HIGHLIGHT_COLOR, atoms)
        self.atom_highlights.extend(atoms)
        if render: self.render()

    def label_atoms(self, atoms=None, render=True, **kwargs):
        kwargs['render'] = False
        if atoms is None:
            atoms = self.mol.atoms
        for atom in atoms:
            self.draw_label(atom.position, atom.name, **kwargs)
        if render: self.render()

    def add_click_callback(self, fn):
        assert callable(fn)
        self._callbacks.add(fn)

    def remove_click_callback(self, fn):
        self._callbacks.remove(fn)

    def handle_click(self, trait_name, old, new):
        # TODO: this should handle more than just atoms
        atom = self.mol.atoms[new['index']]
        if self.selection_group:
            self.selection_group.update_selections(self, {'atoms': [atom]})
        for callback in self._callbacks:
            callback(atom)

    def append_frame(self, positions=None, render=True):
        # override base method - we'll handle frames entirely in python
        # this is copied verbatim from molviz, except for the line noted
        if positions is None:
            positions = self.get_positions()

        positions = self._convert_units(positions)
        try:
            positions = positions.tolist()
        except AttributeError:
            pass

        self.num_frames += 1
        self._frame_positions.append(positions)  # only modification from molviz
        self.show_frame(self.num_frames - 1)
        if render: self.render()

    @utils.doc_inherit
    def show_frame(self, framenum, _fire_event=True, **kwargs):
        # override base method - we'll handle frames using self.set_positions
        # instead of any built-in handlers
        if framenum != self.current_frame:
            self.set_positions(self._frame_positions[framenum])
            self.current_frame = framenum
            if _fire_event and self.selection_group:
                self.selection_group.update_selections(self, {'framenum': framenum})

    def handle_selection_event(self, selection):
        """
        Deals with an external selection event
        :param selection:todo
        :return:
        """
        orb_changed = False
        if 'atoms' in selection:
            self.highlight_atoms(selection['atoms'], render=False)

        if 'framenum' in selection:
            if self.frame_change_callback is not None:
                self.frame_change_callback(selection['framenum'])

            # If user is going to scrub through a trajectory with orbitals,
            # they need to be cached first to avoid locking up the interpreter
            # TODO: parallelize this, also add a "cache" button or something
            if ((self.current_orbital is not None) and
                    (self.current_orbital[1] is not None) and
                    (self.current_orbital not in self._cached_orbitals)):
                self._cache_orb_trajectory(self.current_orbital)

            orb_changed = True
            self.show_frame(selection['framenum'], _fire_event=False, render=False)

        if ('orbname' in selection) and (selection['orbname'] != self.current_orbital):
            orb_changed = True
            self.current_orbital = selection['orbname']

        if ('orbital_isovalue' in selection) and (
                    selection['orbital_isovalue'] != self._orbital_kwargs['isoval']):
            orb_changed = True
            self._orbital_kwargs['isoval'] = selection['orbital_isovalue']

        # redraw the orbital if necessary
        if orb_changed:
            if self.current_orbital is None or self.current_orbital[1] is None:
                self.remove_orbitals()
            else:
                self.draw_orbital(self.current_orbital, render=False, **self._orbital_kwargs)

        self.render()

    def _cache_orb_trajectory(self, orbital):
        print 'Computing wavefunction trajectory for orbital %s ...' % str(orbital),
        sys.stdout.flush()
        for iframe in xrange(self.num_frames):
            self.frame_change_callback(iframe)
            temp = self.get_voldata(orbital,
                                    npts=self._orbital_kwargs.get('npts', None),
                                    _framenum=iframe)
        self._cached_orbitals.add(orbital)
        print 'done.'

    def _convert_units(self, obj):
        try:
            return obj.value_in(self.DISTANCE_UNITS)
        except AttributeError:
            return obj
