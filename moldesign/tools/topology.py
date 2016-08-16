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

ATNUM_VALENCE_CHARGE = {6: {3: -1, 4: 0},
                        7: {2: -1, 3: 0, 4: 1},
                        8: {1: -1, 2: 0, 3: 1},
                        9: {1: 0, 2: 1}}

ION_CHARGE = {1: 1,  # H
              11: 1,  # Na+
              19: 1,  # K+
              12: 2,  # Mg 2+
              20: 2,  # Ca 2+
              9: -1,  # F-
              17: -1,  # Cl-
              35: -1,  # Br-
              53: -1}  # I-


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

    # TODO: scrape formal charge data from the PDB chem comp dictionary

    # cache these values in case we fail to assign charges
    totalcharge = mol.charge
    oldcharges = [atom.formal_charge for atom in mol.atoms]

    for atom in mol.atoms:
        if ignore_nonzero and atom.formal_charge != 0:
            continue

        v = atom.valence
        newcharge = None

        if atom.atnum in ATNUM_VALENCE_CHARGE:
            if v in ATNUM_VALENCE_CHARGE[atom.atnum]:
                newcharge = ATNUM_VALENCE_CHARGE[atom.atnum][v]
            else:
                for oldcharge, a in zip(oldcharges, mol.atoms):
                    a.oldcharge = oldcharge
                mol.charge = totalcharge
                raise UnhandledValenceError(atom)
        elif atom in ION_CHARGE and v == 0:
            newcharge = ION_CHARGE[atom.atnum]

        if newcharge is not None:
            mol.charge += newcharge * u.q_e - atom.formal_charge
            atom.formal_charge = newcharge * u.q_e


@toplevel
def add_missing_data(mol):
    """ Add missing hydrogens, bond orders, and formal charges to a structure (often from the PDB)

    Specifically, this is a convenience function that runs:
    ``mdt.guess_bond_orders``, ``mdt.add_hydrogen``, and ``mdt.assign_formal_charges``

    Note:
        This does NOT add missing residues to biochemical structures. This functionality will be
        available as :meth:`moldesign.add_missing_residues`

    Args:
        mol (moldesign.Molecule): molecule to clean

    Returns:
        moldesign.Molecule: cleaned version of the molecule
    """
    m = mdt.add_hydrogen(mdt.guess_bond_orders(mol))
    assign_formal_charges(m)
    return m


@toplevel
def guess_histidine_states(mol):
    """ Attempt to assign protonation states to histidine residues.

    Note:
        This function is highly unlikely to give accurate results! It is intended for convenience
        when histidine states can easily be guessed from already-present hydrogens or when they are
        judged to be relatively unimportant.

    This can be done simply by renaming HIS residues:
      1. If HE2 and HD1 are present, the residue is renamed to HIP
      2. If only HE2 is present, the residue is renamed to HIE
      3. Otherwise, the residue is renamed to HID (the most common form)

    Args:
        mol (moldesign.Molecule): molecule to change (in place)
    """
    for residue in mol.residues:
        if residue.resname == 'HIS':
            oldname = str(residue)
            if 'HE2' in residue and 'HD1' in residue:
                residue.resname = 'HIP'
            elif 'HE2' in residue:
                residue.resname = 'HIE'
            else:
                residue.resname = 'HID'
            print 'Renaming %s from HIS to %s' % (oldname, residue.resname)







