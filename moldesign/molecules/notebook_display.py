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
import traitlets
import ipywidgets as ipy
import IPython.display

from moldesign import data, utils, uibase, viewer
from moldesign import widgets as mdtwidgets
from moldesign import units as u


class AtomNotebookMixin(object):
    """ Notebook-specific methods for atoms.

    Note:
        This is a mixin class designed only to be mixed into the :class:`Atom` class. Routines
        are separated are here for code organization only - they could be included in the main
        Atom class without changing any functionality
    """
    # @utils.args_from(mdt.molecule.Molecule.draw2d, allexcept=['highlight_atoms'])  # import order
    def draw2d(self, **kwargs):
        """ Draw a 2D viewer with this atom highlighted (Jupyter only).
        In biomolecules, only draws the atom's residue.

        Args:
            width (int): width of viewer in pixels
            height (int): height of viewer in pixels

        Returns:
            mdt.ChemicalGraphViewer: viewer object
        """
        if self.molecule:
            if self.molecule.is_small_molecule:
                return self.molecule.draw2d(highlight_atoms=[self], **kwargs)
            elif self.molecule.is_biomolecule:
                return self.residue.draw2d(highlight_atoms=[self], **kwargs)
            else:
                raise ValueError('No drawing routine specified')
        else:
            raise ValueError('No drawing routine specified')

    # @utils.args_from(mdt.molecule.Molecule.draw2d, allexcept=['highlight_atoms'])  # import order
    def draw3d(self, **kwargs):
        """ Draw a 3D viewer with this atom highlighted (Jupyter only).

        Args:
            width (int): width of viewer in pixels
            height (int): height of viewer in pixels

        Returns:
            mdt.GeometryViewer: viewer object
        """
        return self.molecule.draw3d(highlight_atoms=[self], **kwargs)

    def draw(self, width=300, height=300):
        """ Draw a 2D and 3D viewer with this atom highlighted (notebook only)

        Args:
            width (int): width of viewer in pixels
            height (int): height of viewer in pixels

        Returns:
            ipy.HBox: viewer object
        """
        viz2d = self.draw2d(width=width, height=height, display=False)
        viz3d = self.draw3d(width=width, height=height, display=False)
        return ipy.HBox([viz2d, viz3d])

    def markdown_summary(self):
        """Return a markdown-formatted string describing this atom

        Returns:
            str: markdown-formatted string
        """
        if self.molecule is None:
            lines = ["<h3>Atom %s</h3>" % self.name]
        else:
            lines = ["<h3>Atom %s (index %d)</h3>" % (self.name, self.index)]

        lines.append('**Atomic number**: %d' % self.atnum)
        lines.append("**Mass**: %s" % self.mass)
        lines.append('**Formal charge**: %s' % self.formal_charge)

        if self.molecule is not None:
            lines.append('\n')
            if self.molecule.is_biomolecule:
                if self.pdbindex is not None:
                    lines.append('**PDB serial #**: %s'%self.pdbindex)
                lines.append("**Residue**: %s (index %d)" % (self.residue.name, self.residue.index))
                lines.append("**Chain**: %s" % self.chain.name)
            lines.append("**Molecule**: %s" % self.molecule.name)
            for ibond, (nbr, order) in enumerate(self.bond_graph.iteritems()):
                lines.append('**Bond %d** (order = %d): %s (index %s) in %s' % (
                    ibond + 1, order, nbr.name, nbr.index, nbr.residue.name))

        if self.basis_functions:
            lines.append('**Basis functions:**<br>' + '<br>'.join(map(str,self.basis_functions)))

        if self.ff:
            lines.append('\n**Forcefield partial charge**: %s' % self.ff.partial_charge)
            # TODO: deal with other LJ types, e.g., AB?
            lines.append(u'**Forcefield LJ params**: '
                         u'\u03C3=%s, \u03B5=%s' % (
                             self.ff.ljsigma.defunits(),
                             self.ff.ljepsilon.defunits()))

        # position and momentum
        table = utils.MarkdownTable('', 'x', 'y', 'z')

        table.add_line(['**position /** {}'.format(u.default.length)] +
                       ['%12.3f' % x.defunits_value() for x in self.position])
        table.add_line(['**momentum /** {}'.format(u.default.momentum)] +
                       ['%12.3e' % m.defunits_value() for m in self.momentum])

        if self.molecule is not None and 'forces' in self.molecule.properties:
            table.add_line(['**force /** {.units}'.format(self.force.defunits())] +
                           ['%12.3e' % m.defunits_value() for m in self.force])

        lines.append('\n\n' + table.markdown() + '\n\n')
        # All other assigned properties

        return '<br>'.join(lines)

    def _repr_markdown_(self):
        return self.markdown_summary()


class AtomGroupNotebookMixin(object):
    def draw(self, width=500, height=500, show_2dhydrogens=None, display=False):
        """ Visualize this molecule (Jupyter only).

        Creates a 3D viewer, and, for small molecules, a 2D viewer).

        Args:
            width (int): width of the viewer in pixels
            height (int): height of the viewer in pixels
            show_2dhydrogens (bool): whether to show the hydrogens in 2d (default: True if there
                   are 10 or less heavy atoms, false otherwise)
            display (bool): immediately display this viewer

        Returns:
            moldesign.ui.SelectionGroup
        """

        viz2d = None
        if self.num_atoms < 40:

            viz2d = self.draw2d(width=width, height=height,
                                display=False,
                                show_hydrogens=show_2dhydrogens)
            viz3d = self.draw3d(width=width, height=height,
                                display=False)
            traitlets.link((viz3d, 'selected_atom_indices'), (viz2d, 'selected_atom_indices'))
            views = ipy.HBox([viz2d, viz3d])
        else:
            views = self.draw3d(display=False)

        atom_inspector = uibase.components.AtomInspector()
        traitlets.directional_link(
            (viz2d or views, 'selected_atom_indices'),
            (atom_inspector, 'value'),
            lambda selected_atom_indices: atom_inspector.indices_to_value(selected_atom_indices, self.atoms)
        )

        displayobj = uibase.SelectionGroup([views, atom_inspector])

        if display:
            IPython.display.display(displayobj)
        return displayobj

    def draw3d(self, highlight_atoms=None, **kwargs):
        """ Draw this object in 3D. Jupyter only.

        Args:
            highlight_atoms (List[Atom]): atoms to highlight when the structure is drawn

        Returns:
            mdt.GeometryViewer: 3D viewer object
        """
        self.viz3d = viewer.GeometryViewer(self, **kwargs)
        if highlight_atoms is not None:
            self.viz3d.highlight_atoms(highlight_atoms)
        return self.viz3d

    def draw2d(self, highlight_atoms=None, show_hydrogens=None, **kwargs):
        """
        Draw this object in 2D. Jupyter only.

        Args:
            highlight_atoms (List[Atom]): atoms to highlight when the structure is drawn
            show_hydrogens (bool): whether to draw the hydrogens or not (default: True if there
                   are 10 or less heavy atoms, false otherwise)

        Returns:
            mdt.ChemicalGraphViewer: 2D viewer object
        """
        if show_hydrogens is None:
            show_hydrogens = len(self.heavy_atoms) <= 10
        if not show_hydrogens:
            alist = [atom for atom in self.atoms if atom.atnum > 1]
        else:
            alist = self
        self.viz2d = viewer.DistanceGraphViewer(alist, **kwargs)
        if highlight_atoms: self.viz2d.highlight_atoms(highlight_atoms)
        return self.viz2d


class ResidueNotebookMixin(object):
    """ Mixin class for Residues with notebook-relevenat classes
    """

    def _repr_markdown_(self):
        return self.markdown_summary()

    def markdown_summary(self):
        """ Markdown-formatted information about this residue

        Returns:
            str: markdown-formatted string
        """
        if self.type == 'placeholder':
            return '`%s`' % repr(self)

        if self.molecule is None:
            lines = ["<h3>Residue %s</h3>" % self.name]
        else:
            lines = ["<h3>Residue %s (index %d)</h3>" % (self.name, self.index)]

        if self.type == 'protein':
            lines.append('**Residue codes**: %s / %s' % (self.resname, self.code))
        else:
            lines.append("**Residue code**: %s" % self.resname)
        lines.append('**Type**: %s' % self.type)
        if self.resname in data.RESIDUE_DESCRIPTIONS:
            lines.append('**Description**: %s' % data.RESIDUE_DESCRIPTIONS[self.resname])

        lines.append('**<p>Chain:** %s' % self.chain.name)

        lines.append('**PDB sequence #**: %d' % self.pdbindex)

        terminus = None
        if self.type == 'dna':
            if self.is_3prime_end:
                terminus = "3' end"
            elif self.is_5prime_end:
                terminus = "5' end"
        elif self.type == 'protein':
            if self.is_n_terminal:
                terminus = 'N-terminus'
            elif self.is_c_terminal:
                terminus = 'C-terminus'
        if terminus is not None:
            lines.append('**Terminal residue**: %s of chain %s' % (terminus, self.chain.name))

        if self.molecule is not None:
            lines.append("**Molecule**: %s" % self.molecule.name)

        lines.append("**<p>Number of atoms**: %s" % self.num_atoms)
        if self.backbone:
            lines.append("**Backbone atoms:** %s" % ', '.join(x.name for x in self.backbone))
            lines.append("**Sidechain atoms:** %s" % ', '.join(x.name for x in self.sidechain))
        else:
            lines.append("**Atom:** %s" % ', '.join(x.name for x in self.atoms))

        return '<br>'.join(lines)


class MolNotebookMixin(object):
    """ Methods for displaying molecular information in Jupyter notebooks

    Note:
        This is a mixin class designed only to be mixed into the :class:`Molecule` class. Routines
        are separated are here for code organization only - they could be included in the main
        Molecule class without changing any functionality
    """
    def draw_orbitals(self, **kwargs):
        """ Visualize any calculated molecular orbitals (Jupyter only).

        Returns:
            mdt.orbitals.OrbitalViewer
        """
        if 'wfn' not in self.properties:
            self.calculate_wfn()
        return mdtwidgets.OrbitalViewer(self, **kwargs)


    def markdown_summary(self):
        """A markdown description of this molecule.

        Returns:
            str: Markdown"""
        # TODO: remove leading underscores for descriptor-protected attributes
        lines = ['### Molecule: "%s" (%d atoms)' % (self.name, self.natoms)]

        description = self.metadata.get('description', None)
        if description is not None:
            description = self.metadata.description[:5000]
            url = self.metadata.get('url', None)
            if url is not None:
                description = '<a href="%s" target="_blank">%s</a>' % \
                              (url, description)
            lines.append(description)

        lines.extend([
                 '**Mass**: {:.2f}'.format(self.mass),
                 '**Formula**: %s' % self.get_stoichiometry(html=True),
                 '**Charge**: %s' % self.charge])

        if self.energy_model:
            lines.append('**Potential model**: %s' % str(self.energy_model))

        if self.integrator:
            lines.append('**Integrator**: %s' % str(self.integrator))

        if self.is_biomolecule:
            lines.extend(self.biomol_summary_markdown())

        return '\n\n'.join(lines)

    def _repr_markdown_(self):
        return self.markdown_summary()

    def biomol_summary_markdown(self):
        """A markdown description of biomolecular structure.

        Returns:
            str: Markdown string"""
        lines = []
        if len(self.residues) > 1:
            table = self.get_residue_table()
            lines.append('### Residues')
            # extra '|' here may be workaround for a bug in ipy.markdown?
            lines.append(table.markdown(replace={0: ' '}) + '|')

            lines.append('### Biopolymer chains')
            seqs = []
            for chain in self.chains:
                seq = chain._get_sequence(_html=True)
                if not seq.strip():  # don't write anything if there's no sequence
                    continue

                seqs.append('**%s**: %s' % (chain.name, seq))
            lines.append('<br>'.join(seqs))

        return lines

    def get_residue_table(self):
        """Creates a data table summarizing this molecule's primary structure.

        Returns:
            moldesign.utils.MarkdownTable"""
        table = utils.MarkdownTable(*(['chain'] +
                                      'protein dna rna unknown water solvent'.split()))
        for chain in self.chains:
            counts = {}
            unk = []
            for residue in chain.residues:
                cat = residue.type
                if cat == 'unknown':
                    unk.append(residue.name)
                counts[cat] = counts.get(cat, 0) + 1
            counts['chain'] = '<b>%s</b>' % chain.name
            if 0 < len(unk) <= 4:
                counts['unknown'] = ','.join(unk)
            table.add_line(counts)
        return table

    def get_stoichiometry(self, html=False):
        """ Return this molecule's stoichiometry

        Returns:
            str
        """
        counts = {}
        for atom in self.atoms:
            counts[atom.symbol] = counts.get(atom.symbol, 0) + 1

        my_elements = sorted(counts.keys())
        if html: template = '%s<sub>%d</sub>'
        else: template = '%s%d'
        return ''.join([template % (k, counts[k]) for k in my_elements])


class TrajNotebookMixin(object):
    @utils.kwargs_from(mdtwidgets.trajectory.TrajectoryViewer)
    def draw3d(self, **kwargs):
        """TrajectoryViewer: create a trajectory visualization

        Args:
            **kwargs (dict): kwargs for :class:`moldesign.widgets.trajectory.TrajectoryViewer`
        """
        self._viz = mdtwidgets.trajectory.TrajectoryViewer(self, **kwargs)
        return self._viz
    draw = draw3d  # synonym for backwards compatibility

    def draw_orbitals(self, align=True):
        """ Visualize trajectory with molecular orbitals

        Args:
            align (bool): Align orbital phases (i.e., multiplying by -1 as needed) to prevent sign
               flips between frames

        Returns:
            TrajectoryOrbViewer: create a trajectory visualization
        """
        for frame in self:
            if 'wfn' not in frame:
                raise ValueError("Can't draw orbitals - orbital information missing in at least "
                                 "one frame. It must be calculated with a QM method.")

        if align: self.align_orbital_phases()
        self._viz = mdtwidgets.trajectory.TrajectoryOrbViewer(self)
        return self._viz

    def plot(self, x, y, **kwargs):
        """ Create a matplotlib plot of property x against property y

        Args:
            x,y (str): names of the properties
            **kwargs (dict): kwargs for :meth:`matplotlib.pylab.plot`

        Returns:
            List[matplotlib.lines.Lines2D]: the lines that were plotted

        """
        from matplotlib import pylab
        xl = yl = None
        if type(x) is str:
            strx = x
            x = getattr(self, x)
            xl = '%s / %s' % (strx, x.units)
        if type(y) is str:
            stry = y
            y = getattr(self, y)
            yl = '%s / %s' % (stry, y.units)
        plt = pylab.plot(x, y, **kwargs)
        pylab.xlabel(xl); pylab.ylabel(yl); pylab.grid()
        return plt