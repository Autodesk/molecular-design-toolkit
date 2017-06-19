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
import future.utils

import numpy as np
import imp

import moldesign.units as u
from .. import compute
from ..utils import if_not_none, redirect_stderr
from .. import orbitals
from ..utils import exports

if future.utils.PY2:
    from cStringIO import StringIO
else:
    from io import StringIO


try:
    imp.find_module('pyscf')
except (ImportError, OSError) as exc:
    print('PySCF not installed; using remote docker container')
    force_remote = True
else:
    force_remote = False


@exports
def mol_to_pyscf(mol, basis, symmetry=None, charge=0, positions=None):
    """Convert an MDT molecule to a PySCF "Mole" object"""
    from pyscf import gto
    pyscfmol = gto.Mole()

    positions = if_not_none(positions, mol.positions)
    pyscfmol.atom = [[atom.elem, pos.value_in(u.angstrom)]
                     for atom, pos in zip(mol.atoms, positions)]
    pyscfmol.basis = basis
    pyscfmol.charge = charge
    if symmetry is not None:
        pyscfmol.symmetry = symmetry
    with redirect_stderr(StringIO()) as builderr:
        pyscfmol.build()
    builderr.seek(0)
    for line in builderr:
        if line.strip() == 'Warn: Ipython shell catchs sys.args':
            continue
        else:
            print('PYSCF: ' + line)
    return pyscfmol


# PYSCF appears to have weird names for spherical components?
SPHERICAL_NAMES = {'y^3': (3, -3), 'xyz': (3, -2), 'yz^2': (3, -1), 'z^3': (3, 0),
                   'xz^2': (3, 1), 'zx^2': (3, 2), 'x^3': (3, 3),
                   'x2-y2': (2, 2)}
SPHERICAL_NAMES.update(orbitals.ANGULAR_NAME_TO_COMPONENT)

# TODO: need to handle parameters for max iterations,
# level shifts, requiring convergence, restarts, initial guesses


class StatusLogger(object):
    LEN = 15

    def __init__(self, description, columns, logger):
        self.logger = logger
        self.description = description
        self.columns = columns
        self._init = False
        self._row_format = ("{:<%d}" % self.LEN) + ("{:>%d}" % self.LEN) * (len(columns) - 1)

    def __call__(self, info):
        if not self._init:
            self.logger.status('Starting energy model calculation: %s' % self.description)
            self.logger.status(self._row_format.format(*self.columns))
            self.logger.status(self._row_format.format(*['-' * (self.LEN - 2) for i in self.columns]))
            self._init = True
        self.logger.status(self._row_format.format(*[info.get(c, 'n/a') for c in self.columns]))


@compute.runsremotely(enable=force_remote)
def get_eris_in_basis(basis, orbs):
    """ Get electron repulsion integrals transformed in (in form eri[i,j,k,l] = (ij|kl))
    """
    from pyscf import ao2mo

    pmol = mol_to_pyscf(basis.wfn.mol, basis=basis.basisname)
    eri = ao2mo.full(pmol, orbs.T, compact=True) * u.hartree
    eri.defunits_inplace()
    return orbitals.ERI4FoldTensor(eri, orbs)


@compute.runsremotely(enable=force_remote)
def basis_values(mol, basis, coords, coeffs=None, positions=None):
    """ Calculate the orbital's value at a position in space

    Args:
        mol (moldesign.Molecule): Molecule to attach basis set to
        basis (moldesign.orbitals.BasisSet): set of basis functions
        coords (Array[length]): List of coordinates (with shape ``(len(coords), 3)``)
        coeffs (Vector): List of ao coefficients (optional; if not passed, all basis fn
                 values are returned)

    Returns:
        Array[length]: if ``coeffs`` is not passed, an array of basis fn values at each
               coordinate. Otherwise, a list of orbital values at each coordinate
    """
    from pyscf.dft import numint

    # TODO: more than just create the basis by name ...
    pmol = mol_to_pyscf(mol, basis=basis.basisname, positions=positions)
    aovals = numint.eval_ao(pmol, np.ascontiguousarray(coords.value_in(u.bohr)))
    if coeffs is None:
        return aovals
    else:
        return aovals.dot(coeffs)


