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
# limitations under the License.from __future__ import print_function
import itertools
import string
import collections

import numpy as np

import moldesign as mdt
from . import toplevel
from .. import helpers, utils
from .. import units as u
from ..compute import DummyJob
from ..exceptions import NotCalculatedError
from ..min.base import MinimizerBase
from . import PrimaryStructure, AtomGroup, Bond, HasResidues, BondGraph, MolecularProperties
from ..widgets import WidgetMethod
from .coord_arrays import *


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
        self.constraints.append(mdt.geom.FixedPosition(atom, value=pos))
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
        self.constraints.append(
            mdt.geom.DistanceConstraint(atom1, atom2, value=dist))
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
        self.constraints.append(
            mdt.geom.AngleConstraint(atom1, atom2, atom3, value=angle))
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
        self.constraints.append(
            mdt.geom.DihedralConstraint(atom1, atom2, atom3, atom4, value=angle))
        self._reset_methods()
        return self.constraints[-1]

    def constrain_hbonds(self, usecurrent=False):
        """ Constrains the lengths of all bonds involving hydrogen.

        By default, the lengths will be constrained to their forcefield equilibrium values.
        They can also be constrained to their current values by setting ``usecurrent=True``

        Args:
            mol (moldesign.Molecule): Constrain all h-bonds in this molecule
            usecurrent (bool): if False (default), set the constraint values to their forcefield
               equilibrium values (this will fail if no forcefield is assigned). If True, constrain
               the hydrogen bonds at the current values.

        Raises:
            AttributeError: if usecurrent=False but no forcefield is assigned
        """
        self.constraints.append(mdt.geom.HBondsConstraint(self, usecurrent))
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
        """
        df = self.ndims
        if self.integrator is not None:
            if self.integrator.params.get('remove_translation', False):
                df -= 3
            if self.integrator.params.get('remove_rotation', False):
                if self.num_atoms > 2:
                    df -= 3

        const_hbonds = const_water = False
        if self.integrator is not None:
            const_hbonds = self.integrator.params.get('constrain_hbonds', False)
            const_water = self.integrator.params.get('constrain_water', False)

        if const_hbonds:  # TODO: deal with molecular hydrogen
            for atom in self.atoms:
                if atom.atnum == 1:
                    df -= 1

        if const_water:
            for residue in self.residues:
                if residue.type == 'water':  # constrained water has 6 degrees of freedom
                    if const_hbonds:  # two are already accounted for
                        df -= 1
                    else:
                        df -= 3

        for constraint in self.constraints:  # TODO: deal better with overlaps!
            if (const_hbonds
                    and isinstance(constraint, mdt.geom.DistanceConstraint)
                    and (constraint.a1.atnum == 1 or constraint.a2.atnum == 1)):
                continue
            df -= constraint.dof
        return df

    @property
    def charge(self):
        return self._charge

    @charge.setter
    def charge(self, val):
        if not hasattr(val, 'units'):
            val = val * u.q_e
        self._charge = val

    @property
    def num_electrons(self):
        """int: The number of electrons in the system, based on the atomic numbers and self.charge"""
        return sum(atom.atnum for atom in self.atoms) - self.charge.value_in(u.q_e)

    @property
    def homo(self):
        """int: The array index (0-based) of the highest occupied molecular orbital (HOMO).

        Note:
            This assumes a closed shell ground state! """
        return self.num_electrons // 2 - 1

    @property
    def lumo(self):
        """int: The array index (0-based) of the lowest unoccupied molecular orbital (LUMO).

        Note:
            This assumes a closed shell ground state! """
        return self.num_electrons // 2

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


class MolTopologyMixin(object):
    """ Functions for building and keeping track of bond topology and biochemical structure.

    Note:
        This is a mixin class designed only to be mixed into the :class:`Molecule` class. Routines
        are here for code organization only - they could be included in the main
        Atom class without changing any functionality
    """
    def copy(self, name=None):
        """ Create a copy of the molecule and all of its substructures, metadata, and methods

        Returns:
            Molecule: copied molecule
        """
        if name is None:
            name = self.name + ' copy'
        newmol = Molecule(self.atoms,
                          name=name,
                          pdbname=self.pdbname,
                          charge=self.charge,
                          metadata=self.metadata)
        if self.energy_model is not None:
            newmodel = self._copy_method('energy_model')
            newmol.set_energy_model(newmodel)
        if self.integrator is not None:
            newintegrator = self._copy_method('integrator')
            newmol.set_integrator(newintegrator)
        if self.ff is not None:
            self.ff.copy_to(newmol)
        newmol.constraints = [c.copy(newmol) for c in self.constraints]
        newmol.properties = self.properties.copy_to(mol=newmol)
        return newmol

    def _copy_method(self, methodname):
        method = getattr(self, methodname)
        newmethod = method.__class__()
        newmethod.params.clear()
        newmethod.params.update(method.params)
        return newmethod

    def _rebuild_from_atoms(self):
        """ Rebuild component data structures based on atomic data
        """
        self.is_biomolecule = False
        self.ndims = 3 * self.num_atoms
        self._positions = np.zeros((self.num_atoms, 3)) * u.default.length
        self._momenta = np.zeros((self.num_atoms, 3)) * u.default.momentum
        self.masses = np.zeros(self.num_atoms) * u.default.mass
        self.dim_masses = u.broadcast_to(self.masses, (3, self.num_atoms)).T
        self._assign_atom_indices()
        self.chains.rebuild_hierarchy()
        self._dof = None
        self._topology_changed()

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

    def is_identical(self, other, verbose=False):
        """ Test whether two molecules are "identical"

        We specifically test these quantities for equality:
         - positions
         - momenta
         - chain names
         - residue names
         - atom names / numbers / masses
         - bonds

        Note:
            This tests geometry and topology only; it does not test
            energy models or any calculated properties; it also ignores the ``time`` attribute.

        Args:
            other (moldesign.Molecule): molecule to test against
            verbose (bool): when returning False, print the reasons why


        Returns:
            bool: true if all tested quantities are equal
        """
        if not self.same_topology(other, verbose=verbose):
            return False

        if (self.positions != other.positions).any():
            if verbose:
                mismatch = (self.positions != other.positions).any(axis=1).sum()
                print("%d atoms have different positions." % mismatch)
            return False

        if (self.momenta != other.momenta).any():
            if verbose:
                mismatch = (self.momenta != other.momenta).any(axis=1).sum()
                print("%d atoms have different momenta."%mismatch)
            return False

        return True

    residues = utils.Alias('chains.residues')

    @property
    def num_residues(self):
        return len(self.chains.residues)
    nresidues = numresidues = num_residues

    @property
    def num_chains(self):
        return len(self.chains)
    nchains = numchains = num_chains

    def combine(self, *others):
        """ Create a new molecule from a group of other AtomContainers

        Notes:
            - Chain IDs and sequence numbers are automatically assigned if they are missing
            - Chains will be renamed to prevent chain ID clashes
            - Residue resnames are not changed.

        Args:
            *others (AtomContainer or AtomList or List[moldesign.Atom]):

        Returns:
            mdt.Molecule: a new Molecule that's the union of this structure with all
              others. Chains will be renamed as necessary to avoid clashes.
        """
        new_atoms = []
        charge = 0
        names = []

        chain_names = collections.OrderedDict((x, None) for x in string.ascii_uppercase)
        taken_names = set()
        seen_chains = set()

        for obj in itertools.chain([self], others):
            objatoms = mdt.helpers.get_all_atoms(obj).copy()
            for atom in objatoms:
                chain = atom.chain
                if chain not in seen_chains:
                    seen_chains.add(chain)
                    if chain.pdbindex is None or chain.pdbindex in taken_names:
                        chain.pdbindex = chain.name = next(iter(chain_names.keys()))
                    chain_names.pop(chain.pdbindex, None)
                    taken_names.add(chain.pdbindex)

            new_atoms.extend(objatoms)
            charge += getattr(obj, 'charge', 0*u.q_e)
            if hasattr(obj, 'name'):
                names.append(obj.name)
            elif objatoms[0].molecule is not None:
                names.append('%d atoms from %s' % (len(objatoms), objatoms[0].molecule.name))
            else:
                names.append('list of %d unowned atoms' % len(objatoms))

        return mdt.Molecule(new_atoms,
                            copy_atoms=True,
                            charge=charge,
                            name='%s extended with %d atoms' %
                                 (self.name, len(new_atoms) - self.num_atoms),
                            metadata=utils.DotDict(description=
                                                   'Union of %s' % ', '.join(names)))

    def same_topology(self, other, verbose=False):
        """ Test whether two molecules have equivalent topologies

        We specifically test these quantities for equality:
         - chain names
         - residue names
         - atom names / numbers / masses
         - bonds

        Note:
            This tests geometry and topology only; it does not test
            energy models or any calculated properties; it also ignores the ``time`` attribute.

        Args:
            other (moldesign.Molecule): molecule to test against
            verbose (bool): when returning False, print the reasons why

        Returns:
            bool: true if all tested quantities are equal
        """
        if self.num_atoms != other.num_atoms:
            if verbose: print('INFO: Different numbers of atoms')
            return False

        if self.num_chains != other.num_chains:
            if verbose: print('INFO: Different numbers of chains')
            return False

        if self.num_residues != other.num_residues:
            if verbose: print('INFO: Different numbers of residues')
            return False

        for a1, a2 in itertools.zip_longest(
                itertools.chain(self.atoms, self.residues, self.chains),
                itertools.chain(other.atoms, other.residues, other.chains)):
            if a1.name != a2.name:
                if verbose:
                    print('INFO: %s[%d]: names "%s" and "%s"' % (a1.__class__.__name__, a1.index,
                                                           a1.name, a2.name))
                return False

        if (self.masses != other.masses).any():
            return False

        for a1, a2 in zip(self.atoms, other.atoms):
            if a1.atnum != a2.atnum:
                if verbose:
                    print('INFO: atoms[%d]: atom numbers %d and %d' % (a1.index, a1.atnum, a2.atnum))
                return False

        return self.same_bonds(other, verbose=verbose)

    def same_bonds(self, other, verbose=False):
        for myatom, otheratom in zip(self.atoms, other.atoms):
            mybonds = self.bond_graph[myatom]
            otherbonds = other.bond_graph[otheratom]
            if len(mybonds) != len(otherbonds):
                if verbose:
                    print('INFO: atoms[%d] has %d bonds in self, %d bonds in other' % (
                        myatom.index, len(mybonds), len(otherbonds)))
                return False

            for mynbr, myorder in mybonds.items():
                othernbr = other.atoms[mynbr.index]
                if othernbr not in otherbonds or otherbonds[othernbr] != myorder:
                    if verbose:
                        print('INFO: atoms[%d] bonded to atom[%d] (order %d) in self but not other' % (
                            myatom.index, mynbr.index, myorder))
                    return False
        return True

    def not_identical(self, other):
        return not self.is_identical(other)


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
            raise AttributeError('Cannot simulate; no integrator set for %s' % self)

        init_time = self.time
        traj = self.integrator.run(run_for)
        print('Done - integrated "%s" from %s to %s' % (self, init_time, self.time))
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
            raise AttributeError('Cannot calculate properties; no energy model set for %s' % self)

        if requests is None:
            requests = []

        # Figure out what needs to be calculated,
        # and either launch the job or set the result
        to_calculate = set(list(requests) + self.energy_model.DEFAULT_PROPERTIES)
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

        if 'charge' not in model.params:
            if self.charge != 0:
                model.params.charge = self.charge
        elif model.params['charge'] != self.charge:
            print("Warning: molecular charge (%s) does not match energy model's charge (%s)" % (
                self.charge.value_in(u.q_e), model.params.charge.value_in(u.q_e)))
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
            raise AttributeError('Cannot minimize molecule; no energy model set for %s' % self)

        try:
            trajectory = self.energy_model.minimize(**kwargs)
        except NotImplementedError:
            trajectory = mdt.minimize(self, **kwargs)
        print('Reduced energy from %s to %s' % (trajectory.potential_energy[0],
                                                trajectory.potential_energy[-1]))
        if assert_converged:
            raise NotImplementedError()

        return trajectory

    def _reset_methods(self):
        """
        Called whenever a property is changed that the energy model and/or integrator
        need to know about
        """
        # TODO: what should this do with the property object?
        if self.energy_model is not None:
            self.energy_model._prepped = False
        if self.integrator is not None:
            self.integrator._prepped = False

    def _topology_changed(self):
        self._reset_methods()
        self.ff = None


@toplevel
class Molecule(AtomGroup,
               MolConstraintMixin,
               MolPropertyMixin,
               MolTopologyMixin,
               MolSimulationMixin,
               HasResidues):
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
        copy_atoms (bool): Create the molecule with *copies* of the passed atoms
            (they will be copied automatically if they already belong to another molecule)
        pdbname (str): Name of the PDB file
        charge (units.Scalar[charge]): molecule's formal charge
        electronic_state_index (int): index of the molecule's electronic state
        metadata (dict): Arbitrary metadata dictionary


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
        bond_graph (mdt.molecules.BondGraph): symmetric dictionary specifying bonds between the
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
            - :class:`moldesign.notebook_display.MolNotebookMixin`
            - :class:`moldesign.molecules.MolSimulationMixin`
            - :class:`moldesign.molecules.MolTopologyMixin`
            - :class:`moldesign.molecules.MolConstraintMixin`
            - :class:`moldesign.molecules.MolNotebookMixin`
    """
    positions = ProtectedArray('_positions')
    momenta = ProtectedArray('_momenta')

    draw_orbitals = WidgetMethod('molecules.draw_orbitals')
    _PERSIST_REFERENCES = True  # relevant for `pyccc` RPC calls

    def __init__(self, atomcontainer,
                 name=None,
                 copy_atoms=False,
                 pdbname=None,
                 charge=None,
                 metadata=None):
        super().__init__()

        # initial property init
        self.name = 'uninitialized molecule'
        self._constraints = None
        self._charge = None
        self._properties = None

        atoms, name = self._get_initializing_atoms(atomcontainer, name, copy_atoms)

        if metadata is None:
            metadata = getattr(atomcontainer, 'metadata', utils.DotDict())

        self.atoms = atoms
        self.bond_graph = BondGraph(self)
        self.time = getattr(atomcontainer, 'time', 0.0 * u.default.time)
        self.pdbname = pdbname
        self.constraints = []
        self.energy_model = None
        self.integrator = None
        self.metadata = metadata
        self.electronic_state_index = 0

        if charge is not None:
            self.charge = charge
        else:
            self.charge = getattr(atomcontainer, 'charge',
                                  u.unitsum(atom.formal_charge for atom in self.atoms))

        # Builds the internal memory structures
        self.chains = PrimaryStructure(self)
        self._rebuild_from_atoms()

        if name is not None:
            self.name = name
        elif not self.is_small_molecule:
            self.name = 'unnamed macromolecule'
        else:
            self.name = self.get_stoichiometry()

        self._properties = MolecularProperties(self)
        self.ff = None

    def _get_initializing_atoms(self, atomcontainer, name, copy_atoms):
        """ Make a copy of the passed atoms as necessary, return the name of the molecule
        """
        # copy atoms from another object (i.e., a molecule)
        oldatoms = helpers.get_all_atoms(atomcontainer)
        if copy_atoms or (oldatoms[0].molecule is not None):
            atoms = oldatoms.copy_atoms()
            if name is None:  # Figure out a reasonable name
                if oldatoms[0].molecule is not None:
                    name = oldatoms[0].molecule.name+' copy'
                elif hasattr(atomcontainer, 'name') and isinstance(atomcontainer.name, str):
                    name = utils.if_not_none(name, atomcontainer.name+' copy')
                else:
                    name = 'unnamed'
        else:
            atoms = oldatoms
        return atoms, name

    def __repr__(self):
        try:
            return '<%s (%s), %d atoms>' % (self.name,
                                            self.__class__.__name__,
                                            len(self.atoms))
        except (KeyError, AttributeError):
            return '<molecule (error in __repr__) at %s>' % id(self)

    def __str__(self):
        return 'Molecule: %s' % self.name

    @property
    def constraints(self):
        return self._constraints

    @constraints.setter
    def constraints(self, val):
        self._constraints = utils.ExclusiveList(val, key=utils.methodcaller('_constraintsig'))

    def _repr_markdown_(self):
        """A markdown description of this molecule.

        Returns:
            str: Markdown
        """
        # TODO: remove leading underscores for descriptor-protected attributes
        lines = ['### Molecule: "%s" (%d atoms)'%(self.name, self.natoms)]

        description = self.metadata.get('description', None)
        if description is not None:
            description = self.metadata.description[:5000]
            url = self.metadata.get('url', None)
            if url is not None:
                description = '<a href="%s" target="_blank">%s</a>'% \
                              (url, description)
            lines.append(description)

        lines.extend([
            '**Mass**: %s' % self.mass,
            '**Formula**: %s' % self.get_stoichiometry(html=True),
            '**Charge**: %s' % self.charge])

        if self.energy_model:
            lines.append('**Potential model**: %s'%str(self.energy_model))

        if self.integrator:
            lines.append('**Integrator**: %s'%str(self.integrator))

        lines.extend(self.chains._repr_markdown_())

        return '\n\n'.join(lines)

    def new_bond(self, a1, a2, order):
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
        return sum(atom.nbonds for atom in self.atoms) // 2

    nbonds = num_bonds

    def add_atom(self, newatom):
        """  Add a new atom to the molecule

        Args:
            newatom (moldesign.Atom): The atom to add
                (it will be copied if it already belongs to a molecule)
        """
        self.add_atoms([newatom])

    def add_atoms(self, newatoms):
        """Add new atoms to this molecule.

        *Copies* of the passed atoms will be added if they already belong to another molecule.

        Args:
           newatoms (List[moldesign.Atom]))
        """
        owner = newatoms[0].molecule
        for atom in newatoms:
            if atom.molecule is not owner:
                raise ValueError('Cannot add atoms from multiple sources - add them separately.')

        if owner is not None:  # copy if the atoms are already owned
            newatoms = mdt.AtomList(newatoms).copy()

        self.atoms.extend(newatoms)
        self.bond_graph._add_atoms(newatoms)
        self._rebuild_from_atoms()

    def delete_bond(self, bond_or_atom, a2=None):
        """ Remove this bond from the molecule's topology

        Args:
            bond_or_atom (Bond or Atom): bond to remove (or the first of two atom)
            a2 (Atom): second atom in bond (if first argument was also an atom)
        """
        if a2 is None:  # it's a bond
            a1, a2 = bond_or_atom.a1, bond_or_atom.a2
        else:  # they passed 2 atoms
            a1 = bond_or_atom
        self.bond_graph[a1].pop(a2)
        self._topology_changed()

    def write(self, filename=None, **kwargs):
        """ Write this molecule to a string or file.

        Calls :ref:`moldesign.converters.write`

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
