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

import numpy as np

import moldesign as mdt
from moldesign import helpers, utils
from moldesign.exceptions import NotCalculatedError
from moldesign import units as u
from moldesign.compute import DummyJob
from moldesign.min.base import MinimizerBase

from . import toplevel, Residue, Chain, Instance, AtomContainer, Bond
from .coord_arrays import *


@toplevel
class MolecularProperties(utils.DotDict):
    """ Stores property values for a molecule.
    These objects will be generally created and updated by EnergyModels, not by users.
    """
    def __init__(self, mol, **properties):
        """Initialization: ``properties`` MUST include positions.

        Args:
            mol (Molecule): molecule that these properties are associated with
            **properties (dict): values of molecular properties (MUST include positions as a key)
        """
        # ADD_FEATURE: always return stored properties in the default unit systems
        super(MolecularProperties, self).__init__(positions=mol.positions.copy(), **properties)

    def geometry_matches(self, mol):
        """Returns:
            bool: True if the molecule's ``position`` is the same as these properties' ``position``
        """
        return np.array_equal(self.positions, mol.positions)


class MolConstraintMixin(object):
    """ Functions for applying and managing geometrical constraints.

    Note:
        This is a mixin class designed only to be mixed into the :class:`Molecule` class. Routines
        are separated are here for code organization only - they could be included in the main
        Molecule class without changing any functionality
    """
    def clear_constraints(self):
        """
        Clear all geometry constraints from the molecule.

        Note:
            This does NOT clear integrator options - such as "constrain H bonds"
        """
        self.constraints.clear()
        self._reset_methods()

    def constrain_atom(self, atom, pos=None):
        """  Constrain the position of an atom

        Args:
            atom (moldesign.Atom): The atom to constrain
            pos (moldesign.units.MdtQuantity): position to fix this atom at (default: atom.position) [length]

        Returns:
            moldesign.geometry.FixedPosition: constraint object
        """
        from moldesign import geom
        self.constraints.append(geom.FixedPosition(atom, value=pos))
        self._reset_methods()
        return self.constraints[-1]

    def constrain_distance(self, atom1, atom2, dist=None):
        """  Constrain the distance between two atoms

        Args:
            atom1 (moldesign.Atom)
            atom2 (moldesign.Atom)
            dist ([length]): distance value (default: current distance)

        Returns:
            moldesign.geometry.DistanceConstraint: constraint object
        """
        from moldesign import geom
        self.constraints.append(
            geom.constraints.DistanceConstraint(atom1, atom2, value=dist))
        self._reset_methods()
        return self.constraints[-1]

    def constrain_angle(self, atom1, atom2, atom3, angle=None):
        """  Constrain the bond angle atom1-atom2-atom3

        Args:
            atom1 (moldesign.Atom)
            atom2 (moldesign.Atom)
            atom3 (moldesign.Atom)
            angle ([angle]): angle value (default: current angle)

        Returns:
            moldesign.geometry.AngleConstraint: constraint object
        """
        from moldesign import geom
        self.constraints.append(
            geom.constraints.AngleConstraint(atom1, atom2, atom3, value=angle))
        self._reset_methods()
        return self.constraints[-1]

    def constrain_dihedral(self, atom1, atom2, atom3, atom4, angle=None):
        """  Constrain the bond angle atom1-atom2-atom3

        Args:
            atom1 (moldesign.Atom)
            atom2 (moldesign.Atom)
            atom3 (moldesign.Atom)
            atom4 (moldesign.Atom)
            angle ([angle]): angle value (default: current angle)

        Returns:
            moldesign.geom.AngleConstraint: constraint object
        """
        from moldesign import geom
        self.constraints.append(
            geom.constraints.DihedralConstraint(atom1, atom2, atom3, atom4, value=angle))
        self._reset_methods()
        return self.constraints[-1]


class MolPropertyMixin(object):
    """ Functions for calculating and accessing molecular properties.

    Note:
        This is a mixin class designed only to be mixed into the :class:`Molecule` class. Routines
        are separated are here for code organization only - they could be included in the main
        Molecule class without changing any functionality
    """
    @property
    def mass(self):
        """ u.Scalar[mass]: the molecule's mass
        """
        return sum(self.atoms.mass)

    @property
    def kinetic_energy(self):
        r""" u.Scalar[energy]: Classical kinetic energy :math:`\sum_{\text{atoms}} \frac{p^2}{2m}`
        """
        return helpers.kinetic_energy(self.momenta, self.dim_masses)

    @property
    def kinetic_temperature(self):
        r""" [temperature]: temperature calculated using the equipartition theorem,

        :math:`\frac{2 E_{\text{kin}}}{k_b f}`,

        where :math:`E_{\text{kin}}` is the kinetic energy and :math:`f` is the number of
        degrees of freedom (see :meth:`dynamic_dof <Molecule.dynamic_dof>`)
        """
        return helpers.kinetic_temperature(self.kinetic_energy,
                                           self.dynamic_dof)

    @property
    def dynamic_dof(self):
        """ int: Count the number of spatial degrees of freedom of the system,
        taking into account any constraints

        Note:
            If there are other DOFs not taken into account here, this quantity can be
            set explicitly
        """
        if self._dof is not None:
            return self._dof
        df = self.ndims
        if self.integrator is not None:
            if self.integrator.params.get('remove_translation', False):
                df -= 3
            if self.integrator.params.get('remove_rotation', False):
                if self.num_atoms > 2:
                    df -= 2

        const_hbonds = const_water = False
        if self.integrator is not None:
            const_hbonds = self.integrator.params.get('constrain_hbonds', False)
            const_water = self.integrator.params.get('constrain_water', False)

        if const_hbonds:  # TODO: deal with molecular hydrogen
            for atom in self.atoms:
                if atom.atnum == 1: df -= 1

        if const_water:
            for residue in self.residues:
                if residue.type == 'water':  # constrained water has 6 degrees of freedom
                    if const_hbonds:  # two are already accounted for
                        df -= 1
                    else:
                        df -= 3

        for constraint in self.constraints:  # TODO: deal with more double-counting cases
            if const_hbonds:
                if isinstance(constraint, mdt.geom.DistanceConstraint):
                    # don't double-count constrained hbonds
                    if constraint.a1.atnum == 1 or constraint.a2.atnum == 1: continue
            df -= constraint.dof
        return df

    @dynamic_dof.setter
    def dynamic_dof(self, val):
        self._dof = val

    @property
    def num_electrons(self):
        """int: The number of electrons in the system, based on the atomic numbers and self.charge"""
        return sum(self.atoms.atnum) - self.charge.value_in(u.q_e)

    @property
    def homo(self):
        """int: The array index (0-based) of the highest occupied molecular orbital (HOMO).

        Note:
            This assumes a closed shell ground state! """
        return self.num_electrons/2-1

    @property
    def lumo(self):
        """int: The array index (0-based) of the lowest unoccupied molecular orbital (LUMO).

        Note:
            This assumes a closed shell ground state! """
        return self.num_electrons/2

    @property
    def wfn(self):
        """ moldesign.orbitals.ElectronicWfn: return the molecule's current electronic state,
        if calculated.

        Raises:
            NotCalculatedError: If the electronic state has not yet been calculated at this
                geometry
        """
        return self.get_property('wfn')

    def calc_property(self, name, **kwargs):
        """ Calculate the given property if necessary and return it

        Args:
            name (str): name of the property (e.g. 'potential_energy', 'forces', etc.)

        Returns:
            object: the requested property
        """
        result = self.calculate(requests=[name], **kwargs)
        return result[name]

    def get_property(self, name):
        """ Return the given property if already calculated; raise NotCalculatedError otherwise

        Args:
            name (str): name of the property (e.g. 'potential_energy', 'forces', etc.)

        Raises:
            NotCalculatedError: If the molecular property has not yet been calculated at this
                geometry

        Returns:
            object: the requested property
        """
        if name in self.properties and np.array_equal(self.properties.positions, self.positions):
            return self.properties[name]
        else:
            raise NotCalculatedError(
                    ("The '{0}' property hasn't been calculated yet. "
                     "Calculate it with the molecule.calculate_{0}() method").format(name))

    def calculate_forces(self, **kwargs):
        """ Calculate forces and return them

        Returns:
            units.Vector[force]
        """
        return self.calc_property('forces')

    def calculate_potential_energy(self, **kwargs):
        """ Calculate potential energy and return it

        Returns:
            units.Scalar[energy]: potential energy at this position
        """
        return self.calc_property('potential_energy')

    def calculate_dipole(self, **kwargs):
        """ Calculate forces and return them

        Returns:
            units.Vector[length*charge]: dipole moment at this position (len=3)
        """
        return self.calc_property('dipole')

    def calculate_wfn(self, **kwargs):
        """ Calculate forces and return them

        Returns:
            moldesign.orbitals.ElectronicWfn: electronic wavefunction object
        """
        return self.calc_property('wfn')

    def update_properties(self, properties):
        """
        This is intended mainly as a callback for long-running property calculations.
        When they are finished, they can call this method to update the molecule's properties.

        Args:
            properties (dict): properties-like object. MUST contain a 'positions' attribute.
        """
        if self.properties is None:
            self.properties = properties
        else:
            assert (self.positions == properties.positions).all(), \
                'The molecular geometry does not correspond to these properties'
            self.properties.update()

    @property
    def potential_energy(self):
        """ units.Scalar[energy]: return the molecule's current potential energy, if calculated.

        Raises:
            NotCalculatedError: If the potential energy has not yet been calculated at this
                geometry
        """
        return self.get_property('potential_energy')

    @property
    def forces(self):
        """ units.Vector[force]: return the current force on the molecule, if calculated.

        Raises:
            NotCalculatedError: If the forces have not yet been calculated at this geometry
        """
        return self.get_property('forces')

    @property
    def dipole(self):
        """ units.Vector[length*charge]: return the molecule's dipole moment, if calculated (len=3).

        Raises:
            NotCalculatedError: If the dipole moment has not yet been calculated at this
                geometry
        """
        return self.get_property('dipole')

    @property
    def properties(self):
        """MolecularProperties: Molecular properties calculated at this geometry
        """
        # ADD_FEATURE: some sort of persistent caching so that they aren't lost
        if not self._properties.geometry_matches(self):
            self._properties = MolecularProperties(self)
        return self._properties

    @properties.setter
    def properties(self, val):
        """ Sanity checks - make sure that these properties correspond to the correct geoemtry.
        """
        assert val.geometry_matches(self), \
            "Can't set properties - they're for a different molecular geometry"
        self._properties = val

    # synonyms for backwards compatibility
    calc_wfn = calculate_wfn
    calc_dipole = calculate_dipole
    calc_potential_energy = calculate_potential_energy
    calc_forces = calculate_forces


class MolDrawingMixin(object):
    """ Methods for visualizing molecular structure.

    See Also:
        :class:`moldesign.molecules.atomcollections.AtomContainer`

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
        from moldesign.widgets.orbitals import OrbitalViewer
        if 'wfn' not in self.properties:
            self.calculate_wfn()
        return OrbitalViewer(self, **kwargs)


class MolReprMixin(object):
    """ Methods for creating text-based representations of the molecule

    Note:
        This is a mixin class designed only to be mixed into the :class:`Molecule` class. Routines
        are separated are here for code organization only - they could be included in the main
        Molecule class without changing any functionality
    """
    def __repr__(self):
        try:
            return '<%s (%s), %d atoms>' % (self.name,
                                            self.__class__.__name__,
                                            len(self.atoms))
        except:
            return '<molecule (error in __repr__) at %s>' % id(self)

    def __str__(self):
        return 'Molecule: %s' % self.name

    def markdown_summary(self):
        """A markdown description of this molecule.

        Returns:
            str: Markdown"""
        # TODO: remove leading underscores for descriptor-protected attributes
        lines = ['### Molecule: "%s" (%d atoms)' % (self.name, self.natoms),
                 '**Mass**: {:.2f}'.format(self.mass),
                 '**Formula**: %s' % self.get_stoichiometry(html=True),
                 '**Charge**: %s'%self.charge]

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
                seq = chain.sequence
                if not seq.strip():  # don't write anything if there's no sequence
                    continue

                # deal with extra-long sequences
                seqstring = []
                for i in xrange(0, len(seq), 80):
                    seqstring.append(seq[i:i + 80])
                seqstring = '\n'.join(seqstring)
                seqs.append('**%s**: `%s`' % (chain.name, seqstring))
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
            counts['chain'] = '<pre><b>%s</b></pre>' % chain.name
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
        for symbol in self.atoms.symbol:
            counts[symbol] = counts.get(symbol, 0) + 1

        my_elements = sorted(counts.keys())
        if html: template = '%s<sub>%d</sub>'
        else: template = '%s%d'
        return ''.join([template % (k, counts[k]) for k in my_elements])


class MolTopologyMixin(object):
    """ Functions for building and keeping track of bond topology and biochemical structure.

    Note:
        This is a mixin class designed only to be mixed into the :class:`Molecule` class. Routines
        are separated are here for code organization only - they could be included in the main
        Atom class without changing any functionality
    """
    def copy(self, name=None):
        """ Create a copy of the molecule and all of its substructures

        Returns:
            Molecule: copied molecule

        Note:
            Assigned energy models and integrators are not currently copied, although properties are
        """
        if name is None:
            name = self.name + ' copy'
        newmol = Molecule(self.atoms,
                          name=name,
                          pdbname=self.pdbname,
                          charge=self.charge)
        newmol.properties = self.properties.copy()
        return newmol

    def to_json(self):
        js = mdt.chemjson.jsonify(self,
                                  ('time residues atoms name'
                                   'properties energy_model integrator').split())
        js['chains'] = list(self.chains)
        js['bonds'] = list(self.bonds)
        return js


    def assert_atom(self, atom):
        """If passed an integer, just return self.atoms[atom].
         Otherwise, assert that the atom belongs to this molecule"""
        if type(atom) is int:
            atom = self.mol.atoms[atom]
        else:
            assert atom.molecule is self, "Atom %s does not belong to %s" % (atom, self)
        return atom

    def rebuild(self):
        self.chains = Instance(molecule=self)
        self.residues = []
        self._rebuild_topology()

    def _rebuild_topology(self, bond_graph=None):
        """ Build the molecule's bond graph based on its atoms' bonds

        Args:
            bond_graph (dict): graph to build the bonds from
        """
        if bond_graph is None:
            self.bond_graph = self._build_bonds(self.atoms)
        else:
            self.bond_graph = bond_graph

        self.is_biomolecule = False
        self.ndims = 3 * self.num_atoms
        self._positions = np.zeros((self.num_atoms, 3)) * u.default.length
        self._momenta = np.zeros((self.num_atoms, 3)) * u.default.momentum
        self.masses = np.zeros(self.num_atoms) * u.default.mass
        self.dim_masses = u.broadcast_to(self.masses, (3, self.num_atoms)).T  # TODO: pickling
        self._assign_atom_indices()
        self._assign_residue_indices()
        self._dof = None

    @staticmethod
    def _build_bonds(atoms):
        """ Build a bond graph describing bonds between this list of atoms

        Args:
            atoms (List[moldesign.atoms.Atom])
        """
        # TODO: check atom parents
        bonds = {}

        # First pass - create initial bonds
        for atom in atoms:
            assert atom not in bonds, 'Atom appears twice in this list'
            if hasattr(atom, 'bonds') and atom.bond_graph is not None:
                bonds[atom] = atom.bond_graph
            else:
                bonds[atom] = {}

        # Now make sure both atoms have a record of their bonds
        for atom in atoms:
            for nbr in bonds[atom]:
                if atom in bonds[nbr]:
                    assert bonds[nbr][atom] == bonds[atom][nbr]
                else:
                    bonds[nbr][atom] = bonds[atom][nbr]
        return bonds

    def _assign_atom_indices(self):
        """
        Create geometry-level information based on constituent atoms, and mark the atoms
        as the property of this molecule
        """
        idim = 0
        for idx, atom in enumerate(self.atoms):
            atom._set_molecule(self)
            atom.index = idx
            idim += 3
            self.masses[idx] = atom.mass
            # Here, we index the atom arrays directly into the molecule
            atom._index_into_molecule('_position', self.positions, idx)
            atom._index_into_molecule('_momentum', self.momenta, idx)

    def _assign_residue_indices(self):
        """
        Set up the chain/residue/atom hierarchy
        """
        # TODO: consistency checks

        if self._defchain is None:
            self._defchain = Chain(name=None,
                                   index=None,
                                   molecule=None)

        if self._defres is None:
            self._defres = Residue(name=None,
                                   index=None,
                                   pdbindex=None,
                                   pdbname=None,
                                   chain=self._defchain,
                                   molecule=None)
            self._defchain.add(self._defres)

        default_residue = self._defres
        default_chain = self._defchain
        num_biores = 0

        for atom in self.atoms:
            # if atom has no chain/residue, assign defaults
            if atom.residue is None:
                atom.residue = default_residue
                atom.chain = default_chain
                atom.residue.add(atom)

            # assign the chain to this molecule if necessary
            if atom.chain.molecule is None:
                atom.chain.molecule = self
                atom.chain.index = len(self.chains)

                assert atom.chain.name not in self.chains
                self.chains.add(atom.chain)
            else:
                assert atom.chain.molecule is self

            # assign the residue to this molecule
            if atom.residue.molecule is None:
                atom.residue.molecule = self
                atom.residue.index = len(self.residues)
                self.residues.append(atom.residue)
                if atom.residue.type in ('dna', 'rna', 'protein'): num_biores += 1
            else:
                assert atom.chain.molecule is self

        self.is_biomolecule = (num_biores >= 2)
        self.nchains = self.n_chains = self.num_chains = len(self.chains)
        self.nresidues = self.n_residues = self.num_residues = len(self.residues)


class MolSimulationMixin(object):
    """ Functions calculating energies, running dynamics, and minimizing geometry.

    Note:
        This is a mixin class designed only to be mixed into the :class:`Molecule` class. Routines
        are separated are here for code organization only - they could be included in the main
        Atom class without changing any functionality
    """
    def run(self, run_for):
        """ Starts the integrator's default integration

        Args:
            run_for (int or [time]): number of steps or amount of time to run for

        Returns:
            moldesign.trajectory.Trajectory
        """
        if self.integrator is None:
            raise ValueError('Cannot simulate; no integrator set for %s' % self)

        init_time = self.time
        traj = self.integrator.run(run_for)
        print 'Done - integrated "%s" from %s to %s' % (self, init_time, self.time)
        return traj

    def calculate(self, requests=None, wait=True, use_cache=True):
        """
        Runs a potential energy calculation on the current geometry, returning the requested quantities.
        If `requests` is not passed, the properties specified in the energy_models DEFAULT_PROPERTIES
        will be calculated.

        Args:
            requests (List[str]): list of quantities for the model to calculate,
                e.g. ['dipole', 'forces']
            wait (bool): if True, wait for the calculation to complete before returning. \
                 If false, return a job object - this will not update the molecule's properties!
            use_cache (bool): Return cached results if possible

        Returns:
            MolecularProperties
        """
        if self.energy_model is None:
            raise ValueError('Cannot calculate properties; no energy model set for %s' % self)

        if requests is None:
            requests = []

        # Figure out what needs to be calculated,
        # and either launch the job or set the result
        to_calculate = set(requests + self.energy_model.DEFAULT_PROPERTIES)
        if use_cache:
            to_calculate = to_calculate.difference(self.properties)
        if len(to_calculate) == 0:
            job = self.properties
        else:
            job = self.energy_model.calculate(to_calculate)

        if wait:
            # We'll wait for the job to complete, then
            # returns the molecule's calculated properties
            if hasattr(job, 'wait'):
                job.wait()
                properties = job.result
            else:
                properties = job
            self.properties.update(properties)
            return self.properties
        else:
            # We're not waiting for the job to complete - return a job object
            if hasattr(job, 'wait'):
                return job
            else:
                return DummyJob(job)

    def set_energy_model(self, model, **params):
        """ Associate an energy model with this molecule

        Args:
            model (moldesign.methods.EnergyModelBase): The energy model to associate with this
                molecule
            **params (dict): a dictionary of parameters for the model

        Note:
            For convenience, ``model`` can be an instance, a class, or a constructor (with
            call signature ``model(**params) -> model instance)``
        """
        # If passed the class or constructor, create an instance of the energy model
        if not issubclass(type(model), mdt.models.base.EnergyModelBase) and callable(model):
            model = model()

        self.energy_model = model
        self.properties = MolecularProperties(self)
        model.mol = self
        model.params.update(params)

        if 'charge' in model.params:
            if model.params.charge is None:
               model.params.charge = self.charge
            elif model.params.charge != self.charge:
                print "Warning: molecular charge (%d) does not match energy model's charge (%d)" % (
                    self.charge, model.params.charge)
        model._prepped = False

    def set_integrator(self, integrator, **params):
        """ Associate an integrator with this molecule

        Args:
            integrator (moldesign.integrators.IntegratorBase): The integrator to associate with this
                molecule
            **params (dict): a dictionary of parameters for the integrator

        Note:
            For convenience, ``integrator`` can be an instance, a class, or a constructor (with
            call signature ``integrator(**params) -> integrator instance)``
        """
        # If passed the class or constructor, create an instance of the integrator
        if (not issubclass(type(integrator), mdt.integrators.base.IntegratorBase) and
                callable(integrator)):
            integrator = integrator()

        self.integrator = integrator
        integrator.mol = self
        integrator.params.update(params)
        integrator._prepped = False

    @utils.args_from(MinimizerBase,
                     allexcept=['self'],
                     inject_kwargs={'assert_converged': False})
    def minimize(self, assert_converged=False, **kwargs):
        """Perform an energy minimization (aka geometry optimization or relaxation).

        If ``force_tolerance`` is not specified, the program defaults are used.
        If specified, the largest force component must be less than force_tolerance
        and the RMSD must be less than 1/3 of it. (based on GAMESS OPTTOL keyword)

        Args:
            assert_converged (bool): Raise an exception if the minimization does not converged.

        Returns:
            moldesign.trajectory.Trajectory
        """
        if self.energy_model is None:
            raise ValueError('Cannot minimize molecule; no energy model set for %s' % self)

        try:
            trajectory = self.energy_model.minimize(**kwargs)
        except NotImplementedError:
            trajectory = mdt.minimize(self, **kwargs)
        print 'Reduced energy from %s to %s' % (trajectory.potential_energy[0],
                                                trajectory.potential_energy[-1])
        if assert_converged:
            raise NotImplementedError()

        return trajectory

    def _reset_methods(self):
        """
        Called whenever a property is changed that the energy model and/or integrator
        need to know about
        """
        # TODO: what should this do with the property object?
        # TODO: handle duplicate constraints (this happens a lot, and is bad)
        if self.energy_model is not None:
            self.energy_model._prepped = False
        if self.integrator is not None:
            self.integrator._prepped = False

    def configure_methods(self):
        """ Interactively configure this molecule's simulation methods (notebooks only)

        Returns:
            ipywidgets.Box: configuration widget
        """
        import ipywidgets as ipy

        children = []
        if self.energy_model:
            children.append(self.energy_model.configure())

        if self.integrator:
            children.append(self.integrator.configure())

        return ipy.VBox(children)


@toplevel
class Molecule(AtomContainer,
               MolConstraintMixin,
               MolPropertyMixin,
               MolDrawingMixin,
               MolReprMixin,
               MolTopologyMixin, MolSimulationMixin):
    """
    ``Molecule`` objects store a molecular system, including atoms, 3D coordinates, molecular
    properties, biomolecular entities, and other model-specific information. Interfaces with
    simulation models take place through the molecule object.

    Molecule objects will generally be created by reading files or parsing other input; see, for
    example: :meth:`moldesign.read`, :meth:`moldesign.from_smiles`,
    :meth:`moldesign.from_pdb`, etc.

    This constructor is useful, however for copying other molecular structures (see
    examples below).

    Args:
        atomcontainer (AtomContainer or AtomList or List[moldesign.Atom]): atoms that make up
            this molecule.

            Note:
                If the passed atoms don't already belong to a molecule, they will be assigned
                to this one. If they DO already belong to a molecule, they will be copied,
                leaving the original molecule untouched.

        name (str): name of the molecule (automatically generated if not provided)
        bond_graph (dict): dictionary specifying bonds between the atoms - of the form
            ``{atom1:{atom2:bond_order, atom3:bond_order}, atom2:...}``
            This structure must be symmetric; we require
            ``bond_graph[atom1][atom2] == bond_graph[atom2][atom1]``
        copy_atoms (bool): Create the molecule with *copies* of the passed atoms
            (they will be copied automatically if they already belong to another molecule)
        pdbname (str): Name of the PDB file
        charge (units.Scalar[charge]): molecule's formal charge
        electronic_state_index (int): index of the molecule's electronic state


    Examples:
        Use the ``Molecule`` class to create copies of other molecules and substructures thereof:
        >>> benzene = mdt.from_name('benzene')
        >>> benzene_copy = mdt.Molecule(benzene, name='benzene copy')

        >>> protein = mdt.from_pdb('3AID')
        >>> carbon_copies = mdt.Molecule([atom for atom in protein.atoms if atom.atnum==6])
        >>> first_residue_copy = mdt.Molecule(protein.residues[0])

    **Molecule instance attributes:**

    Attributes:
        atoms (AtomList): List of all atoms in this molecule.
        bond_graph (dict): symmetric dictionary specifying bonds between the
           atoms:

               ``bond_graph = {atom1:{atom2:bond_order, atom3:bond_order}, atom2:...}``

               ``bond_graph[atom1][atom2] == bond_graph[atom2][atom1]``
        residues (List[moldesign.Residue]): flat list of all biomolecular residues in this molecule
        chains (Dict[moldesign.Chain]): Biomolecular chains - individual chains can be
            accessed as ``mol.chains[list_index]`` or ``mol.chains[chain_name]``
        name (str): A descriptive name for molecule
        charge (units.Scalar[charge]): molecule's formal charge
        constraints (List[moldesign.geom.GeometryConstraint]): list of constraints
        ndims (int): length of the positions, momenta, and forces arrays (usually 3*self.num_atoms)
        num_atoms (int): number of atoms (synonym: natoms)
        num_bonds (int): number of bonds (synonym: nbonds)
        positions (units.Array[length]): Nx3 array of atomic positions
        momenta (units.Array[momentum]): Nx3 array of atomic momenta
        masses (units.Vector[mass]): vector of atomic masses
        dim_masses (units.Array[mass]): Nx3 array of atomic masses (for numerical convenience -
           allows you to calculate velocity, for instance, as
           ``velocity = mol.momenta/mol.dim_masses``
        time (units.Scalar[time]): current time in dynamics
        energy_model (moldesign.models.base.EnergyModelBase): Object that calculates
            molecular properties - driven by `mol.calculate()`
        integrator (moldesign.integrators.base.IntegratorBase): Object that drives movement of 3D
            coordinates in time, driven by mol.run()
        is_biomolecule (bool): True if this molecule contains at least 2 biochemical residues


    **Molecule methods and properties**

    See also methods offered by the mixin superclasses:

            - :class:`moldesign.molecules.AtomContainer`
            - :class:`moldesign.molecules.MolPropertyMixin`
            - :class:`moldesign.molecules.MolDrawingMixin`
            - :class:`moldesign.molecules.MolSimulationMixin`
            - :class:`moldesign.molecules.MolTopologyMixin`
            - :class:`moldesign.molecules.MolConstraintMixin`
            - :class:`moldesign.molecules.MolReprMixin`
    """
    positions = ProtectedArray('_positions')
    momenta = ProtectedArray('_momenta')

    def __init__(self, atomcontainer,
                 name=None, bond_graph=None,
                 copy_atoms=False,
                 pdbname=None,
                 charge=None,
                 electronic_state_index=0):
        super(Molecule, self).__init__()

        # copy atoms from another object (i.e., a molecule)
        oldatoms = helpers.get_all_atoms(atomcontainer)

        if copy_atoms or (oldatoms[0].molecule is not None):
            #print 'INFO: Copying atoms into new molecule'
            atoms = oldatoms.copy()
            if name is None:  # Figure out a reasonable name
                if oldatoms[0].molecule is not None:
                    name = oldatoms[0].molecule.name + ' copy'
                elif hasattr(atomcontainer, 'name') and isinstance(atomcontainer.name, str):
                    name = utils.if_not_none(name, atomcontainer.name + ' copy')
                else:
                    name = 'unnamed'
        else:
            atoms = oldatoms

        self.atoms = atoms
        self.time = 0.0 * u.default.time
        self.name = 'uninitialized molecule'
        self._defres = None
        self._defchain = None
        self.pdbname = pdbname
        self.constraints = utils.ExclusiveList(key=utils.methodcaller('_constraintsig'))
        self.energy_model = None
        self.integrator = None
        self.electronic_state_index = electronic_state_index

        if charge is not None:
            self.charge = charge
            if not hasattr(charge, 'units'):  # assume fundamental charge units if not explicitly set
                self.charge *= u.q_e
        else:
            self.charge = sum(atom.formal_charge for atom in self.atoms)

        # Builds the internal memory structures
        self.chains = Instance(molecule=self)
        self.residues = []
        self._rebuild_topology(bond_graph=bond_graph)

        if name is not None:
            self.name = name
        elif not self.is_small_molecule:
            self.name = 'unnamed macromolecule'
        else:
            self.name = self.get_stoichiometry()

        self._properties = MolecularProperties(self)
        self.ff = utils.DotDict()

    def newbond(self, a1, a2, order):
        """ Create a new bond

        Args:
            a1 (moldesign.Atom): First atom in the bond
            a2 (moldesign.Atom): Second atom in the bond
            order (int): order of the bond

        Returns:
            moldesign.Bond
        """
        assert a1.molecule == a2.molecule == self
        return a1.bond_to(a2, order)

    @property
    def velocities(self):
        """ u.Vector[length/time]: Nx3 array of atomic velocities
        """
        return (self.momenta/self.dim_masses).defunits()

    @velocities.setter
    def velocities(self, value):
        self.momenta = value * self.dim_masses

    @property
    def num_bonds(self):
        """int: number of chemical bonds in this molecule"""
        return sum(atom.nbonds for atom in self.atoms)/2

    nbonds = num_bonds

    def addatom(self, newatom):
        """  Add a new atom to the molecule

        Args:
            newatom (moldesign.Atom): The atom to add
                (it will be copied if it already belongs to a molecule)
        """
        self.addatoms([newatom])

    def addatoms(self, newatoms):
        """Add new atoms to this molecule.
        For now, we really just rebuild the entire molecule in place.

        Args:
           newatoms (List[moldesign.Atom]))
        """
        self._reset_methods()

        for atom in newatoms: assert atom.molecule is None
        self.atoms.extend(newatoms)

        # symmetrize bonds between the new atoms and the pre-existing molecule
        bonds = self._build_bonds(self.atoms)
        for newatom in newatoms:
            for nbr in bonds[newatom]:
                if nbr in self.bond_graph:  # i.e., it's part of the original molecule
                    bonds[nbr][newatom] = bonds[newatom][nbr]

        self._rebuild_topology(bonds)

    def deletebond(self, bond):
        """ Remove this bond from the molecule's topology

        Args:
            Bond: bond to remove
        """
        self.bond_graph[bond.a1].pop(bond.a2)
        self.bond_graph[bond.a2].pop(bond.a1)

    def _force_converged(self, tolerance):
        """ Return True if the forces on this molecule:
        1) Are less than tolerance in every dimension
        2) have an RMS of less than 1/3 the tolerance value

        Args:
            tolerance (units.Scalar[force]): force tolerance

        Returns:
            bool: True if RMSD force is less than this quantity
        """
        forces = self.calc_forces()
        if forces.max() > tolerance: return False
        rmsd2 = forces.dot(forces) / self.ndims
        if rmsd2 > tolerance * tolerance / 3.0: return False
        return True

    def write(self, filename=None, **kwargs):
        """ Write this molecule to a string or file.

        This is a convenience method for :ref:`moldesign.converters.write`

        Args:
            filename (str): filename to write (if not passed, write to string)
            format (str): file format (if filename is not passed, format must be specified)
                Guessed from file extension if not passed
        """
        # TODO: make it easier to do the right thing, which is write to .pkl.bz2
        return mdt.write(self, filename=filename, **kwargs)

    @property
    def is_small_molecule(self):
        """bool: True if molecule's mass is less than 500 Daltons (not mutually exclusive with
        :meth:`self.is_biomolecule <Molecule.is_biomolecule>`)"""
        return self.mass <= 500.0 * u.amu

    @property
    def bonds(self):
        """ Iterator over all bonds in the molecule

        Yields:
            moldesign.atoms.Bond: bond object
        """
        for atom in self.bond_graph:
            for nbr in self.bond_graph[atom]:
                if atom.index > nbr.index: continue  # don't double count
                yield Bond(atom, nbr)
