
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
import copy

import moldesign as mdt
from moldesign import units as u


class ForcefieldParams(object):
    """ Stores the forcefield parameters for a specific molecule.

    The current implementation is heavily based on a ParmEd object,
    with fairly thin MDT wrappers providing an interface and handling unit conversions.

    Args:
        mol (mdt.Molecule): Molecule this force field is for
        ffobj (parmed.Structure OR mdt.AmberParms OR Forcefield): forcefield specification

    Attributes:
        parmed_obj (parmed.Structure): parmed object containing the forcefield parameters
    """

    def __init__(self, mol, ffobj):
        import parmed

        self.mol = mol

        if isinstance(ffobj, parmed.Structure):
            self.parmed_obj = copy.copy(ffobj)
            self.sourcedata = ffobj
        elif hasattr(ffobj, 'to_parmed'):
            self.sourcedata = ffobj
            self.parmed_obj = ffobj.to_parmed()
        else:
            raise ValueError('Unrecognized force field class "%s"' % ffobj.__class__.__name__)

    def copy_to(self, mol):
        mol.ff = self.__class__(mol, copy.copy(self.parmed_obj))
        return mol.ff

    def __setstate__(self, state):
        """ Workaround for https://github.com/ParmEd/ParmEd/issues/874
        This function can be removed once parmed 0.7.4 is released - AMV 5.19.17
        """
        self.__dict__.update(state)
        self.parmed_obj.initialize_topology()

    def get_atom_terms(self, atom):
        return AtomTerms(atom, self.parmed_obj.atoms[atom.index])

    def get_bond_term(self, bond_or_atom1, atom2=None):
        if atom2 is None:
            bond = bond_or_atom1
        else:
            bond = mdt.Bond(bond_or_atom1,
                            atom2)

        pmdatom = self.parmed_obj.atoms[bond.a1.index]

        for pmdbond, pmdpartner in zip(pmdatom.bonds, pmdatom.bond_partners):
            if pmdpartner.idx == bond.a2.index:
                return BondTerm(bond,
                                pmdbond)
        else:
            raise ValueError("No ForceField term found for bond: %s" % bond)

    def get_term(self, *atoms):
        pmdterm = self._get_pmd_term(*atoms)
        return LEN_TO_TERM[len(atoms)](atoms, pmdterm)

    def _get_pmd_term(self, *atoms):
        termlist = getattr(self.parmed_obj.atoms[atoms[0].index],
                           LEN_TO_TERM[len(atoms)].pmdlist)
        searchatoms = [self.parmed_obj.atoms[atom.index] for atom in atoms]

        for term in termlist:
            if term.same_atoms(searchatoms):
                return term

        else:
            raise ValueError("No ForceField term found with atoms: %s" % atoms)


class ParmedAtomAttribute(object):
    def __init__(self, attrname, units):
        self.attrname = attrname
        self.units = units

    @staticmethod
    def _get_pmd_obj(instance):
        return instance.pmdobj

    def __get__(self, instance, owner):
        obj = self._get_pmd_obj(instance)
        return getattr(obj, self.attrname) * self.units

    def __set__(self, instance, value):
        obj = self._get_pmd_obj(instance)
        setattr(obj, self.attrname, value.value_in(self.units))


class ParmedTermAttribute(ParmedAtomAttribute):
    @staticmethod
    def _get_pmd_obj(instance):
        return instance.pmdobj.type

    def __set__(self, *args):
        # these are defined by *TYPE*, so will need to create singleton classes
        # when modifying these terms.
        raise NotImplementedError()


class ForceFieldTerm(object):
    """ MDT-like proxy for a forcefield term
    """

    def __init__(self, mdtobj, pmdobj):
        self.mdtobj = mdtobj
        self.pmdobj = pmdobj

    def __str__(self):
        return "%s for %s, ParmEd params %s" % (self.__class__.__name__,
                                                self.mdtobj,
                                                str(self.pmdobj.type)[1:-1])  # strips brackets

    def __repr__(self):
        return '<%s>' % str(self)


class AtomTerms(ForceFieldTerm):
    partial_charge = ParmedAtomAttribute('charge', u.q_e)
    ljsigma = ParmedAtomAttribute('sigma', u.angstrom)
    ljepsilon = ParmedAtomAttribute('epsilon', u.kcalpermol)


class BondTerm(ForceFieldTerm):
    pmdlist = 'bonds'
    force_constant = ParmedTermAttribute('k', u.kcalpermol / u.angstrom)
    equilibrium_length = ParmedTermAttribute('req', u.angstrom)


class AngleTerm(ForceFieldTerm):
    pmdlist = 'angles'
    force_constant = ParmedTermAttribute('k', u.kcalpermol / u.degrees)
    equilibrium_length = ParmedTermAttribute('req', u.angstrom)


class DihedralTerm(ForceFieldTerm):
    pmdlist = 'dihedrals'
    force_constant = ParmedTermAttribute('phi_k', u.kcalpermol / u.degrees)
    periodicity = ParmedTermAttribute('per', 1)
    phase = ParmedTermAttribute('phase', u.degrees)
    coulomb_14_factor = ParmedTermAttribute('scee', 1.0)
    lj_14_factor = ParmedTermAttribute('scnb', 1.0)

LEN_TO_TERM = {1: AtomTerms, 2: BondTerm, 3: AngleTerm, 4: DihedralTerm}


class FFParameters(object):
    """ DEPRACATED: will be removed in favor of ParmEd interface

    This object contains assigned force field parameters for a specific system
    The base focuses on the AMBER / CHARMM - type force
    """
    # TODO: this needs to describe attenuation for things close together
    # TODO: deal with nonbonded exceptions
    def __init__(self, bonds, angles, dihedrals, partial_charges, lennard_jones):
        self.bonds = bonds
        self.angles = angles
        self.dihedrals = dihedrals

        self.partial_charges = partial_charges  # maps atoms to their partial charges
        self.lennard_jones = lennard_jones  # maps atoms to LJ terms

        # lookups
        self.bond_term = {term.bond: term for term in self.bonds}
        self.angle_term = {tuple(term.atoms): term for term in self.angles}
        for term in self.angles: self.angle_term[tuple(reversed(term.atoms))] = term
        self.dihedral_term = {}
        for term in self.dihedrals:
            self.dihedral_term.setdefault(tuple(term.atoms), []).append(term)
            self.dihedral_term.setdefault(tuple(reversed(term.atoms)), []).append(term)

