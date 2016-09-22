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
import itertools
import numpy as np
import moldesign as mdt
from moldesign import units as u

from .constraints import FixedCoordinate, FixedPosition

# TODO: create dynamics wrapper that uses timestep to explicitly calculate constraint forces


def shake_positions(mol, prev_positions, max_cycles=100, use_masses=True):
    """ Satisfy all molecular constraints using the SHAKE algorithm

    Args:
        mol (moldesign.Molecule): molecule with unsatisfied constraints
        prev_positions (mdt.units.Array[length]): positions prior to move
        max_cycles (int): halt and raise an exception after this many cycles

    Note:
        This algorithm does badly if constraint function gradients go to 0.

    References:
        R. Elber, A. Ruymgaart, B. Hess: SHAKE Parallelization.
            Eur Phys J Spec Top. 2011 Nov 1; 200(1): 211.
            doi:10.1140/epjst/e2011-01525-9
    """
    constraints = []
    for c in mol.constraints:  # Replace FixedPosition with 3 FixedCoordinates - it's better behaved
        if isinstance(c, FixedPosition):
            for i in xrange(3):
                vec = np.zeros(3)
                vec[i] = 1.0
                constraints.append(FixedCoordinate(c.atom,  vec, value=c.value[i]))
        else:
            constraints.append(c)

    # ---array shapes---
    # prevgrad, currgrad: (num_constr, 3N)
    # values: (num_constr,)
    # A: (num_constr, num_constr)
    # multipliers: (num_constr, )
    # delta: (3N,)
    curr = mol.positions.copy()
    prev = prev_positions.copy()
    mol.positions = prev
    prevgrad = np.array([_clean_grad_array(c.gradient()) for c in constraints])
    mol.positions = curr

    if use_masses:
        dim_masses = mol.dim_masses
    else:
        dim_masses = np.ones((mol.num_atoms, 3)) * u.default.mass
    flat_masses = dim_masses.defunits_value().flatten()

    for i in xrange(max_cycles):
        for c in mol.constraints:
            if not c.satisfied():
                break
        else:
            return  # e.g., we're done

        # Get constraint derivatives
        # Note: we remove units here because pint does not handle arrays with heterogeneous units
        values = np.array([c.error().defunits_value() for c in constraints])
        currgrad = np.array([_clean_grad_array(c.gradient()) for c in constraints])

        A = np.dot(currgrad/flat_masses, prevgrad.T)
        multipliers = np.linalg.solve(A, values)

        # reapply units and adjust positions
        delta = multipliers.dot(prevgrad).reshape(mol.num_atoms, 3) * (
            u.default.mass * u.default.length)

        mol.positions -= delta/dim_masses

    else:
        raise mdt.ConvergenceFailure('SHAKE did not converge after %d iterations'%
                                     max_cycles)


def _clean_grad_array(a):
    """ Remove units and flatten array
    """
    return a.defunits_value().flatten()








