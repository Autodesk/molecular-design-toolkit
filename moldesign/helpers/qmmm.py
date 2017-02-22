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
import moldesign as mdt

LINKBONDRATIO = 0.709  # fixed ratio of C-C to C-H bond length for link atoms


def create_link_atoms(mol, qmatoms):
    """ Create hydrogen caps for bonds between QM and MM regions.

    Each link atom will have ``metadata.mmatom``, ``metadata.mmpartner`` attributes to identify the
    atom it replaces and the atom it's bonded to in the MM system.

    Raises:
        ValueError: if any MM/QM atom is bonded to more than one QM/MM atom, or the bond
                    order is not one

    Returns:
        List[mdt.Atom]: list of link atoms
    """
    linkatoms = []
    qmset = set(qmatoms)
    for qmatom in qmatoms:
        mmatom = _get_mm_nbr(mol, qmatom, qmset)
        if mmatom is None:
            continue

        la = mdt.Atom(atnum=1, name='HL%d' % len(linkatoms),
                      metadata={'mmatom': mmatom, 'mmpartner': qmatom})
        linkatoms.append(la)

    set_link_atom_positions(linkatoms)
    return linkatoms


def _get_mm_nbr(mol, qmatom, qmset):
    mm_nbrs = [nbr for nbr in qmatom.bonded_atoms
               if nbr not in qmset]
    if len(mm_nbrs) == 0:
        return None

    # everything below is sanity checks
    mmatom = mm_nbrs[0]
    if len(mm_nbrs) != 1:
        raise ValueError('QM atom %s is bonded to more than one MM atom' % qmatom)
    if mol.bond_graph[qmatom][mmatom] != 1:
        raise ValueError('Bond crossing QM/MM boundary (%s - %s) does not have order 1'
                         % (qmatom, mmatom))

    if qmatom.atnum != 6 or mmatom.atnum != 6:
        print ('WARNING: QM/MM bond involving non-carbon atoms: %s - %s' %
               (qmatom, mmatom))
    mm_qm_nbrs = [qmnbr for qmnbr in mmatom.bonded_atoms
                  if qmnbr in qmset]
    if len(mm_qm_nbrs) != 1:
        raise ValueError('MM atom %s is bonded to more than one QM atom'%mmatom)
    return mmatom


def set_link_atom_positions(linkatoms):
    """
    Set link atom positions using a fixed ratio of MM bond length to QM bond length

    Warnings:
        - This is only valid for

        - Presumably, the most "correct" way to do this is to place the hydrogen in order to
          match the force exterted on the QM atom by the MM atom. This is not currently supported.

    Args:
        linkatoms (List[mdt.Atom]): list of link atoms to set positions for

    References:
        http://www.nwchem-sw.org/index.php/Qmmm_link_atoms
    """
    for atom in linkatoms:
        nbr = atom.metadata.mmpartner
        proxy = atom.metadata.mmatom
        dist = LINKBONDRATIO * nbr.distance(proxy)
        atom.position = (nbr.position +
                         dist * mdt.mathutils.normalized(proxy.position - nbr.position))
