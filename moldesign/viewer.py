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
import sys
from itertools import product

import numpy as np
import IPython.display as dsp
import webcolors

from nbmolviz.drivers3d import MolViz_3DMol
from nbmolviz.widget2d import MolViz2DBaseWidget

import moldesign as mdt
import moldesign.units as u
from moldesign import utils
from moldesign.data import color_rotation

from moldesign.geometry import perpendicular


def is_color(s):
    """ Do our best to determine if "s" is a color spec that can be converted to hex
    :param s:
    :return:
    """
    def in_range(i): return 0 <= i <= int('0xFFFFFF', 0)

    try:
        if type(s) == int:
            return in_range(s)
        elif type(s) not in (str, unicode):
            return False
        elif s in webcolors.css3_names_to_hex:
            return True
        elif s[0] == '#':
            return in_range(int('0x' + s[1:], 0))
        elif s[0:2] == '0x':
            return in_range(int(s, 0))
        elif len(s) == 6:
            return in_range(int('0x' + s, 0))
    except ValueError:
        return False


DEF_CATEGORICAL = 'Paired'
DEF_SEQUENTIAL = None  # should be inferno, but that's only MPL >1.5

def colormap(cats, mplmap='auto'):
    # should make it easy to choose one for:
    #  categorical data
    #  sequential (low, high important)
    #  diverging data (low, mid, high important)
    # Can deal with numerical and categorical data
    # we'll treat ints as categories for now
    global DEF_SEQUENTIAL
    from matplotlib import cm

    if hasattr(cm, 'inferno'):
        DEF_SEQUENTIAL = 'inferno'
    else:
        DEF_SEQUENTIAL = 'BrBG'

    # strip units
    units = None
    if hasattr(cats[0], 'magnitude'):
        arr = u.to_units_array(cats)
        units = arr.units
        cats = arr.magnitude

    if not isinstance(cats, np.ndarray) and not isinstance(cats[0], float):  # treat as
        # categorical
        values = np.zeros(len(cats), dtype='float')
        to_int = collections.OrderedDict()
        for i, item in enumerate(cats):
            if item not in to_int:
                to_int[item] = len(to_int)
            values[i] = to_int[item]
        if mplmap == 'auto':
            mplmap = DEF_CATEGORICAL
    else:  # it's numerical
        values = np.array(cats, dtype='float')
        if mplmap == 'auto':
            mplmap = DEF_SEQUENTIAL

    cmap = getattr(cm, mplmap)
    mx = values.max()
    mn = values.min()
    r = (values - mn) / (mx - mn)  # rescale to [0.0,1.0]
    rgb = cmap(r)
    hexcolors = [webcolors.rgb_to_hex(np.array(r[:3]) * 256) for r in rgb]
    return hexcolors


# Right now hard-coded to use 3DMol driver - will need to be configurable when there's more than one
# In theory, the only tie to the specific driver should be the class that we inherit from
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

    def calc_orb_grid(self, orbname, npts=30, padding=2.5 * u.angstrom):
        """
        Calculate orbitals
        :param orbname: Either orbital index (for canonical orbitals) or a tuple ( [orbital type], [orbital index] ) \
         where [orbital type] is a keyword (e.g., canonical, natural, nbo, ao, etc)
        :param npts: number of points in each dimension of the grid
        :param padding: extent of the grid beyond the furthest atom in each dimension
        """
        # TODO: why is orbname sometimes a np.int64 instead of an int?
        # NEWFEATURE: limit bounding box based on the non-zero atomic centers. Useful for localized orbitals,
        #             which otherwise require high resolution
        try:
            orbtype, orbidx = orbname
        except TypeError:
            orbtype = 'canonical'
            orbidx = orbname

        orbs = self.wfn.orbitals[orbtype]  # this is the wfn for the current frame
        return orbs.calculate_orb_grid(orbidx, npoints=npts, padding=padding)

    def get_orbnames(self):
        """
        :return:
        """
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
        # override base method - we'll handle frames using self.set_positions instead of any built-in handlers
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
            if (self.current_orbital is not None) and (
                self.current_orbital not in self._cached_orbitals):
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

        if self.current_orbital is not None and orb_changed:
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


def draw_structure(mol, **kwargs):
    atoms = [atom for atom in mol.atoms if atom.atnum != 1]
    names = [atom.symbol for atom in atoms]
    return ChemicalGraphViewer(atoms, carbon_labels=False, names=names,
                               **kwargs)


class ChemicalGraphViewer(MolViz2DBaseWidget):
    MAXATOMS = 200

    def __init__(self, mol,
                 carbon_labels=True,
                 names=None,
                 display=False,
                 _forcebig=False,
                 **kwargs):

        self.carbon_labels = carbon_labels
        try:
            self.atoms = mol.atoms
        except AttributeError:
            self.atoms = mol
        else:
            self.mol = mol

        if not _forcebig and len(self.atoms) > self.MAXATOMS:
            raise ValueError('Refusing to draw more than 200 atoms in 2D visualization. '
                             'Override this with _forcebig=True')

        if names is None: self.names = [atom.name for atom in self.atoms]

        self.atom_indices = {atom: i for i, atom in enumerate(self.atoms)}
        self.selection_group = None
        self.selection_id = None
        super(ChemicalGraphViewer, self).__init__(self.atoms, **kwargs)
        self.set_click_callback(callback=self.handle_click)
        if display: dsp.display(self)

    def __reduce__(self):
        """These don't get passed around,
        so we send NOTHING"""
        return utils.make_none, tuple()

    def to_graph(self, atoms):
        nodes, links = [], []
        for i1, atom1 in enumerate(atoms):
            nodes.append(dict(atom=atom1.name, index=i1))
            if atom1.atnum == 6 and not self.carbon_labels:
                nodes[-1].update({'atom': '',
                                  'size': 0.5,
                                  'color': 'darkgray'})
            for neighbor, order in atom1.bond_graph.iteritems():
                if neighbor not in self.atom_indices: continue
                nbr_idx = self.atom_indices[neighbor]
                if nbr_idx < i1:
                    links.append({'source': i1,
                                  'target': nbr_idx,
                                  'bond': order})
        graph = dict(nodes=nodes, links=links)
        return graph

    def get_atom_index(self, atom):
        return self.atom_indices[atom]

    def handle_click(self, trait_name, old, new):
        clicked_atoms = [self.atoms[new]]
        if self.selection_group:
            self.selection_group.update_selections(self, {'atoms': clicked_atoms})

    def handle_selection_event(self, selection):
        """
        Deals with an external selection event
        :param selection:
        :return:
        """
        if 'atoms' in selection:
            self.highlight_atoms(
                [a for a in selection['atoms'] if a in self.atom_indices])


class DistanceGraphViewer(ChemicalGraphViewer):
    def __init__(self, atoms,
                 distance_sensitivity=(3.0 * u.ang, 7.0 * u.ang),
                 bond_strength=1.0,
                 angstrom_to_px=22.0,
                 minimum_bond_strength=0.2,
                 nonbond_strength=0.66,
                 charge=-300,
                 **kwargs):
        dmin, dmax = distance_sensitivity
        self.minimum_bond_strength = minimum_bond_strength
        self.dmin = dmin.value_in(u.angstrom)
        self.dmax = dmax.value_in(u.angstrom)
        self.drange = self.dmax - self.dmin
        self.bond_strength = bond_strength
        self.angstrom_to_px = angstrom_to_px
        self.nonbond_strength = nonbond_strength
        self.colored_residues = {}
        kwargs['charge'] = charge
        super(DistanceGraphViewer, self).__init__(atoms,**kwargs)

    def to_graph(self, atoms):
        graph = super(DistanceGraphViewer, self).to_graph(atoms)

        # Deal with covalent bonds
        for link in graph['links']:
            a1 = atoms[link['source']]
            a2 = atoms[link['target']]
            link['strength'] = self.bond_strength
            link['distance'] = a1.distance(a2).value_in(
                u.angstrom) * self.angstrom_to_px

        # Add distance restraints for non-bonded atoms
        for i1, atom1 in enumerate(atoms):
            for i2 in xrange(i1 + 1, len(atoms)):
                atom2 = atoms[i2]
                if atom1 in atom2.bond_graph: continue

                distance = atom1.distance(atom2).value_in(u.angstrom)
                if distance > self.dmax: continue

                strength = self.nonbond_strength * min(
                    float((1.0 - (distance - self.dmin) / self.drange) ** 2),
                    1.0)

                if strength < self.minimum_bond_strength: continue

                link = {'distance': float(distance * self.angstrom_to_px),
                        'source': i1, 'target': i2,
                        'strength': strength, 'bond': 0}
                graph['links'].append(link)

        return graph

    def color_by_residue(self, black_residue=None):
        # Color by residue
        residues = {}
        for atom in self.atoms:
            if atom.residue in residues:
                residues[atom.residue].append(atom)
            else:
                residues[atom.residue] = [atom]
        if black_residue:
            black_atoms = residues.pop(black_residue)
            self.set_atom_style(outline_color='black', atoms=black_atoms)
        for ires, res in enumerate(residues):
            self.set_atom_style(fill_color=color_rotation[ires], atoms=residues[res])
            self.colored_residues[res] = color_rotation[ires]

    def draw_contacts(self, group1, group2, radius=2.25 * u.angstrom,
                      label=True):
        for atom1, atom2 in product(group1, group2):
            if atom1.index == atom2.index: continue
            if atom1 in atom2.bond_graph: continue
            skip = False
            for nbr in atom2.bond_graph:
                if atom1 in nbr.bond_graph:
                    skip = True
            if skip: continue
            dst = atom1.distance(atom2)
            if dst <= radius:
                self.set_bond_style([[atom1, atom2]],
                                    width=1, dash_length=1, opacity=1.0, color='black')
                if label:
                    self.set_bond_label([atom1, atom2],
                                        text='%.1f ang'%dst.value_in('angstrom'),size=8)


def make_contact_view(entity, view_radius=5.0*u.ang,
                      contact_radius=2.25*u.ang,
                      angstrom_to_px=44.0,
                      **kwargs):
    """

    :type entity: moldesign.biounits.Entity
    :param kwargs:
    :return:
    """
    from moldesign.atoms import AtomList

    try:
        focus_atoms = AtomList(entity.atoms)
    except AttributeError:
        focus_atoms = AtomList(entity)

    # get the complete set of atoms to display
    view_atoms = focus_atoms.atoms_within(view_radius)
    all_atoms = AtomList(view_atoms + focus_atoms)
    assert len(set(all_atoms)) == len(view_atoms) + len(focus_atoms)
    viewer = DistanceGraphViewer(all_atoms, **kwargs)
    viewer.color_by_residue(black_residue=focus_atoms[0].residue)
    viewer.draw_contacts(focus_atoms, view_atoms, radius=contact_radius)
    return viewer


class BondSelectorBase(GeometryViewer):
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
        super(BondSelectorBase, self).__init__(mol=mol, render=False, **kwargs)
        self.atom_callbacks = []
        self.bond_callbacks = []
        self.click_callbacks = []
        self.vdw(radius=self.ATOMRADIUS, render=False)
        self.draw_all_bonds(render=True)

    def set_positions(self, *args, **kwargs):
        render = kwargs.get('render',True)
        kwargs['render'] = False
        super(BondSelectorBase, self).set_positions(*args, **kwargs)
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
