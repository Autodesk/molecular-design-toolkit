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
from .. import units as u
from ..mathutils import Eigenspace
from .. import utils


@utils.exports
class PrincipalMomentsOfInertia(Eigenspace):
    """ Calculates the moments of inertia of a molecule at a given position

    Notes:
        This object is NOT updated if the initializing molecule's positions change.
        This allows it to be used to reorient other molecules into the same coordinate system

    Args:
        mol (moldesign.Molecule): molecule to calculate PMIs for
        mass_centered (bool): calculate PMIs relative to molecule's center of mass (default: True)

    Attributes:
        axes (List[Vector[len=3]]]): principal rotation axes
        moments (Vector[length**2 * mass, len=3]): rotational inertias for the three axes
    """
    def __init__(self, mol, mass_centered=True):
        from scipy.linalg import eig
        inertia_tensor, com = get_inertia_tensor(mol, mass_centered=mass_centered)

        evals, evecs = eig(inertia_tensor.defunits_value())  # strip units before going to scipy
        evals = evals*u.default.mass*u.default.length ** 2  # reapply units to evals

        assert (abs(evecs.imag) < 1e-12).all()  # hopefully everything is real?
        assert (abs(evals.imag) < 1e-12 * evals.units).all()

        super().__init__(evals.real, evecs.T)
        self._com = com
        self.sort()

    @property
    def moments(self):
        return self.evals

    @property
    def axes(self):
        return self.evecs

    def reorient(self, mol):
        """ Rotates molecular coordinates into the PMI coordinate system

        The molecule will be translated (if necessary), then rotated so that x lies along PMI 1,
        y along PMI 2, and z along PMI 3.

        Args:
            mol (moldesign.Molecule): molecule to re-orient - its current positinos will be changed
        """
        if self._com is not None:
            mol.com -= self._com
        mol.positions = self.transform(mol.positions)


@utils.exports
def get_inertia_tensor(mol, mass_centered=True):
    """ Calculates the moment of inertia for a molecule

    Args:
        mol (moldesign.Molecule): calculate the inertia tensor of this  molecule
        mass_centered (bool): if True, calcualte tensor relative to molecule's center of mass. If
           False, calculate it relative to the current origin.

    Returns:
        Matrix[mass*length*length]: moment of inertia tensor
    """
    inertia_tensor = np.zeros((3,3)) * u.default.mass * u.default.length**2
    com = mol.com if mass_centered else None

    pos = mol.positions
    if mass_centered:
        pos = pos - mol.com

    for i in range(3):
        inertia_tensor[i,i] = (mol.masses * (pos[:,i-1]**2 + pos[:,i-2]**2)).sum()
        for j in range(i+1,3):
            di = (mol.masses * pos[:,i] * pos[:,j]).sum()
            inertia_tensor[i,j] = -di
            inertia_tensor[j,i] = -di
    return inertia_tensor, com