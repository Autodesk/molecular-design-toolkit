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

from . import MolecularOrbitals
from ..utils import DotDict

class ElectronicWfn(object):
    """ Stores the results of a quantum chemistry calculation.

    This is necessarily pretty flexible, but generally stores an LCAO wfn and one or more sets of
    orbitals. Can also store CI vectors, etc.

    These objects will usually be created by quantum chemical energy models.

    Args:
        mol (moldesign.Molecule): Molecule this wavefunction belongs to
        num_electrons (int): number of electrons in this wavefunction
        model (moldesign.models.base.EnergyModelBase): The model this wavefunction was created with
        aobasis (moldesign.orbitals.BasisSet): The basis functions for the enclosed orbitals
        nbasis (int): number of AO basis functions
        fock_ao (moldesign.units.Array[energy]): fock matrix in the AO basis
        positions (moldesign.units.Array[length]): positions of the nuclei for this wfn
        civectors (np.ndarray): CI vectors (if applicable)
        description (str): text describing the wfn (e.g. 'RHF/STO-3G', 'CAS(2,2)/SA3/6-31G**')
        density_matrix_ao (np.ndarray): density matrix in the ao basis
    """
    def __init__(self, mol, num_electrons,
                 model=None,
                 aobasis=None, fock_ao=None,
                 positions=None,
                 civectors=None,
                 description=None,
                 density_matrix_ao=None):
        self.mol = mol
        self.model = model
        self.civectors = civectors
        self.aobasis = aobasis
        if aobasis:
            self.nbasis = len(self.aobasis)
        else:
            self.nbasis = None
        self.orbitals = DotDict()
        self.fock_ao = fock_ao
        self.num_electrons = num_electrons
        self.homo = self.num_electrons//2 - 1
        self.lumo = self.homo + 1
        self._has_canonical = False
        self.density_matrix_ao = density_matrix_ao
        self.description = description

        if positions is None:
            self.positions = mol.positions.copy()
        else:
            self.positions = positions.copy()

        if self.aobasis is not None:
            self.orbitals['atomic'] = self.aobasis
            self.aobasis.wfn = self
            for orb in self.aobasis.orbitals:
                orb.wfn = self

    def __repr__(self):
        return '<ElectronicWfn (%s) of %s>' % (self.description, str(self.mol))

    def __str__(self):
        return '%s wfn' % self.description

    def set_canonical_mos(self, orbs):
        if orbs.wfn is None:
            orbs.wfn = self

        self._has_canonical = True

    def align_orbital_phases(self, other, assert_same=True):
        """Align this wavefunction's orbitals to have the same phase as those in `other`.

        Args:
            other (ElectronicWfn): reference wfn to align phases with
            assert_same (bool): raise an exception if the two wavefunctions do not have the same
                 kinds of orbitals
        """
        for orbtype in self.orbitals:
            if orbtype not in other.orbitals:
                if assert_same: assert False, '%s has orbital type %s, but %s does not.' % (self, orbtype, other)
                else: continue
            self.orbitals[orbtype].align_phases(other.orbitals[orbtype])

    def run_nbo(self, **kwargs):
        from moldesign.interfaces import nbo_interface
        nbo_interface.run_nbo(self.mol, **kwargs)

    def add_orbitals(self, orbs, orbtype='canonical', **kwargs):
        mo_object = MolecularOrbitals(orbs,
                                      wfn=self,
                                      orbtype=orbtype)
        self.orbitals[orbtype] = mo_object
        if orbtype == 'canonical' and not self._has_canonical:
            self.set_canonical_mos(mo_object)
        return mo_object

    @property
    def molecular_orbitals(self):
        """A synonym for self.orbitals['canonical'], since this is usually what's wanted"""
        return self.orbitals['canonical']

    @molecular_orbitals.setter
    def molecular_orbitals(self, val):
        """A synonym for self.orbitals['canonical'], since this is usually what's wanted"""
        self.orbitals['canonical'] = val
