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

from itertools import product
from IPython import display as dsp

from nbmolviz.widget2d import MolViz2DBaseWidget

import moldesign as mdt
from moldesign import utils
import moldesign.units as u

from . import toplevel, ColorMixin


@toplevel
class ChemicalGraphViewer(MolViz2DBaseWidget, ColorMixin):
    """ Create a JSON-format graph representing the chemical structure and draw it using the
    NBMolViz 2D widget.

    Args:
        mol (moldesign.molecules.AtomContainer): A collection of atoms (eg a list of atoms,
            a residue, a molecule. etc)
        carbon_labels (bool): If True, draw atom names for carbons
        names (List[str]): (optional) a list of strings to label the atoms in the drawing
            (default: ``[atom.name for atom in mol.atoms]``)
        display (bool): immediately display this drawing
    """

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
            self.atoms = mdt.AtomList(mol)
            self.mol = self.atoms
        else:
            self.mol = mol

        if not _forcebig and len(self.atoms) > self.MAXATOMS:
            raise ValueError('Refusing to draw more than 200 atoms in 2D visualization. '
                             'Override this with _forcebig=True')

        if names is None:
            names = []
            for atom in self.atoms:
                if atom.formal_charge == 0:
                    names.append(atom.name)
                else:
                    names.append(atom.name + _charge_str(atom.formal_charge))

        self.names = names

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
            nodes.append(dict(atom=self.names[i1], index=i1))
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
        """ Return the atom's index in this object's storage
        """
        return self.atom_indices[atom]

    def unset_color(self, atoms=None, render=None):
        self.set_color('white', atoms)

    def handle_click(self, trait_name, old, new):
        clicked_atoms = [self.atoms[new]]
        if self.selection_group:
            self.selection_group.update_selections(self, {'atoms': clicked_atoms})

    def handle_selection_event(self, selection):
        """ Highlight atoms in response to a selection event

        Args:
            selection (dict): Selection event from :mod:`moldesign.uibase.selectors`
        """
        if 'atoms' in selection:
            self.highlight_atoms(
                [a for a in selection['atoms'] if a in self.atom_indices])


def _charge_str(q):
    q = q.value_in(u.q_e)
    if q == 0:
        return ''
    elif q == 1:
        return '+'
    elif q == -1:
        return '-'
    elif q > 0:
        return '+%d' % q
    else:
        return str(q)


@toplevel
class DistanceGraphViewer(ChemicalGraphViewer):
    """ Create a 2D graph that includes edges with 3D information. This gives a 2D chemical that
    shows contacts from 3D space.

    Args:
        mol (moldesign.molecules.AtomContainer): A collection of atoms (eg a list of atoms,
            a residue, a molecule. etc)
        distance_sensitivity (Tuple[u.Scalar[length]]): a tuple containing the minimum and
            maximum 3D distances to create edges for (default: ``(3.0*u.ang, 7.0*u.ang)``)
        bond_edge_weight (float): edge weight for covalent bonds
        nonbond_weight_factor (float): scale non-covalent edge weights by this factor
        angstrom_to_px (int): number of pixels per angstrom
        charge (int): the force-directed layout repulsive "charge"
    """
    def __init__(self, atoms,
                 distance_sensitivity=(3.0 * u.ang, 7.0 * u.ang),
                 bond_edge_weight=1.0,
                 minimum_edge_weight=0.2,
                 nonbond_weight_factor=0.66,
                 angstrom_to_px=22.0,
                 charge=-300,
                 **kwargs):
        dmin, dmax = distance_sensitivity
        self.minimum_bond_strength = minimum_edge_weight
        self.dmin = dmin.value_in(u.angstrom)
        self.dmax = dmax.value_in(u.angstrom)
        self.drange = self.dmax - self.dmin
        self.bond_strength = bond_edge_weight
        self.angstrom_to_px = angstrom_to_px
        self.nonbond_strength = nonbond_weight_factor
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
    from moldesign.molecules import AtomList

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