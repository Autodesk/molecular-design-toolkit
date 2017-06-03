"""
This module contains various utility functions that are exposed to API users
"""
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

import moldesign as mdt
from moldesign import units as u

from . import toplevel, __all__ as _pkgall

from moldesign.interfaces.openbabel import add_hydrogen, guess_bond_orders, set_protonation
from moldesign.interfaces.pdbfixer_interface import mutate, add_water
from moldesign.interfaces.ambertools import create_ff_parameters
from moldesign.interfaces.ambertools import calc_am1_bcc_charges, calc_gasteiger_charges

_pkgall.extend(('add_hydrogen guess_bond_orders mutate add_water'
                ' create_ff_parameters calc_am1_bcc_charges calc_gasteiger_charges '
                'set_protonation').split())

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
def set_hybridization_and_ph(mol, ph=7.4):
    """ Add missing hydrogens, bond orders, and formal charges

    Specifically, this is a convenience function that runs:
    ``mdt.guess_bond_orders``, ``mdt.add_hydrogen``, and ``mdt.assign_formal_charges``

    Note:
        This does NOT add missing residues to biochemical structures. This functionality will be
        available as :meth:`moldesign.add_missing_residues`

    Args:
        mol (moldesign.Molecule): molecule to clean
        ph (float): assigned pH. Assign protonation states using the default OpenBabel pKa model

    Returns:
        moldesign.Molecule: cleaned version of the molecule
    """
    m1 = mdt.guess_bond_orders(mol)
    m2 = mdt.add_hydrogen(m1)
    m3 = mdt.set_protonation(m2, ph)
    return m3


@toplevel
def set_hybridization_and_saturate(mol):
    """ Assign bond orders, saturate with hydrogens, and assign formal charges

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
    m1 = mdt.guess_bond_orders(mol)
    m2 = mdt.add_hydrogen(m1)
    assign_formal_charges(m2)
    return m2


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
            print('Renaming %s from HIS to %s' % (oldname, residue.resname))


@toplevel
def split_chains(mol, distance_threshold=1.75*u.angstrom):
    """ Split a molecule's chains into unbroken biopolymers and groups of non-polymers

    This function is non-destructive - the passed molecule will not be modified.

    Specifically, this function will:
       - Split any chain with non-contiguous biopolymeric pieces into single, contiguous polymers
       - Remove any solvent molecules from a chain into their own chain
       - Isolate ligands from each chain into their own chains

    Args:
        mol (mdt.Molecule): Input molecule
        distance_threshold (u.Scalar[length]): if not ``None``, the maximum distance between
           adjacent residues for which we consider them "contiguous". For PDB data, values greater
           than 1.4 Angstrom are eminently reasonable; the default threshold of 1.75 Angstrom is
           purposefully set to be extremely cautious (and still much lower than the distance to
           the *next* nearest neighbor, generally around 2.5 Angstrom)

    Returns:
        mdt.Molecule: molecule with separated chains
    """

    tempmol = mol.copy()

    def bonded(r1, r2):
        if r2 not in r1.bonded_residues:
            return False
        if distance_threshold is not None and r1.distance(r2) > distance_threshold:
            return False
        return True

    def addto(chain, res):
        res.chain = None
        chain.add(res)
        for atom in res:
            atom.chain = chain

    allchains = [mdt.Chain(tempmol.chains[0].name)]
    for chain in tempmol.chains:
        chaintype = chain.residues[0].type
        solventchain = mdt.Chain(None)
        ligandchain = mdt.Chain(None)

        for ires, residue in enumerate(chain.residues):
            if residue.type == 'unknown':
                thischain = ligandchain
            elif residue.type in ('water', 'solvent', 'ion'):
                thischain = solventchain
            else:
                assert residue.type == chaintype
                if ires != 0 and not bonded(residue.prev_residue, residue):
                    allchains.append(mdt.Chain(None))
                thischain = allchains[-1]

            addto(thischain, residue)

        for c in (solventchain, ligandchain):
            if c.num_atoms > 0:
                allchains.append(c)

    return mdt.Molecule(allchains)








