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

"""
This module contains various utility functions that are exposed to API users
"""
import moldesign as mdt
from moldesign import units as u

from . import toplevel, __all__ as _pkgall

from moldesign.interfaces.openbabel import add_hydrogen, guess_bond_orders
from moldesign.interfaces.pdbfixer_interface import mutate, add_water
from moldesign.interfaces.ambertools import assign_forcefield, parameterize
from moldesign.interfaces.ambertools import calc_am1_bcc_charges, calc_gasteiger_charges

_pkgall.extend(('add_hydrogen guess_bond_orders mutate add_water'
                ' assign_forcefield parameterize calc_am1_bcc_charges calc_gasteiger_charges').split())


@toplevel
def assign_formal_charges(mol, ignore_nonzero=True):
    """ Assign formal charges to C,N,O,F atoms in this molecule based on valence

    Args:
        mol (moldesign.Molecule): Molecule to assign formal charges to. The formal charges of its
           atoms and its total charge will be adjusted in place.
        ignore_nonzero (bool): If formal charge is already set to a nonzero value, ignore this atom

    Note:
        This method ONLY applies to C,N, O and F, based on a simple valence model.
        These results should be manually inspected for consistency.

    Raises:
        UnhandledValenceError: for cases not handled by the simple valence model

    References:
        These assignments are illustrated by the formal charge patterns in
            http://www.chem.ucla.edu/~harding/tutorials/formalcharge.pdf
    """
    from moldesign.exceptions import UnhandledValenceError

    # cache these values in case we fail to assign charges
    totalcharge = mol.charge
    oldcharges = [atom.formal_charge for atom in mol.atoms]

    def fail():  # reset all charges before raising exception
        for oldcharge, a in zip(oldcharges, mol.atoms):
            a.oldcharge = oldcharge
        mol.charge = totalcharge
        raise UnhandledValenceError(atom)

    for atom in mol.atoms:
        if ignore_nonzero and atom.formal_charge != 0: continue

        v = atom.valence
        newcharge = None

        if atom.atnum == 6:
            if v == 3:
                newcharge = -1
            elif v == 4:
                newcharge = 0
            else:
                fail()

        elif atom.atnum == 7:
            if v == 2:
                newcharge = -1
            elif v == 3:
                newcharge = 0
            elif v == 4:
                newcharge = 1
            else:
                fail()

        elif atom.atnum == 8:
            if v == 1:
                newcharge = -1
            elif v == 2:
                newcharge = 0
            elif v == 3:
                newcharge = 1
            else:
                fail()

        elif atom.atnum == 9:
            if v == 1:
                newcharge = 0
            elif v == 2:
                newcharge = 1
            else:
                fail()

        if newcharge is not None:
            mol.charge += newcharge * u.q_e - atom.formal_charge
            atom.formal_charge = newcharge * u.q_e

@toplevel
def clean_pdb(mol):
    """ Attempt to clean up a molecule from PDB format that may be missing data

    Specifically, this is a convenience function that runs:
    ``mdt.guess_bond_orders``, ``mdt.add_hydrogen``, and ``mdt.assign_formal_charges``

    Args:
        mol (moldesign.Molecule): molecule to clean

    Returns:
        moldesign.Molecule: cleaned version of the molecule
    """
    m = mdt.add_hydrogen(mdt.guess_bond_orders(mol))
    assign_formal_charges(m)
    return m









