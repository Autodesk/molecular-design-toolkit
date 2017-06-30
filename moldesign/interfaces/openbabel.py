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
import sys
import string

import moldesign.molecules.atomcollections

try:
    import pybel as pb
    import openbabel as ob
    # WARNING: this is the real library, not our interface - this works because of absolute
    # imports. We should probably rename this interface

except ImportError:
    force_remote = True
    sys.stderr.write('Info: OpenBabel not installed; will run in docker container\n')
else:
    force_remote = False

import moldesign as mdt
from ..compute.runsremotely import runsremotely
from .. import units as u
from ..utils import exports


def read_stream(filelike, format, name=None):
    """ Read a molecule from a file-like object

    Note:
        Currently only reads the first conformation in a file

    Args:
        filelike: a file-like object to read a file from
        format (str): File format: pdb, sdf, mol2, bbll, etc.
        name (str): name to assign to molecule

    Returns:
        moldesign.Molecule: parsed result
    """
    molstring = str(filelike.read())  # openbabel chokes on unicode
    return read_string(molstring, format, name=name)


@runsremotely(enable=force_remote)
def read_string(molstring, format, name=None):
    """ Read a molecule from a file-like object

    Note:
        Currently only reads the first conformation in a file

    Args:
        molstring (str): string containing file contents
        format (str): File format: pdb, sdf, mol2, bbll, etc.
        name (str): name to assign to molecule

    Returns:
        moldesign.Molecule: parsed result
    """
    pbmol = pb.readstring(format, molstring)
    mol = pybel_to_mol(pbmol, name=name)
    return mol


@runsremotely(enable=force_remote)
def write_string(mol, format):
    """ Create a file from the passed molecule

    Args:
        mol (moldesign.Molecule): molecule to write
        format (str): File format: pdb, sdf, mol2, bbll, etc.

    Returns:
        str: contents of the file

    References:
        https://openbabel.org/docs/dev/FileFormats/Overview.html
    """
    pbmol = mol_to_pybel(mol)
    if format == 'smi':  # TODO: always kekulize, never aromatic
        outstr = pbmol.write(format=format).strip()
    else:
        outstr = pbmol.write(format=format)
    return str(outstr)


@runsremotely(enable=force_remote)
def guess_bond_orders(mol):
    """Use OpenBabel to guess bond orders using geometry and functional group templates.

    Args:
        mol (moldesign.Molecule): Molecule to perceive the bonds of

    Returns:
        moldesign.Molecule: New molecule with assigned bonds
    """
    # TODO: pH, formal charges
    pbmol = mol_to_pybel(mol)
    pbmol.OBMol.PerceiveBondOrders()
    newmol = pybel_to_mol(pbmol)
    return newmol


@runsremotely(enable=force_remote)
def add_hydrogen(mol):
    """Add hydrogens to saturate atomic valences.

    Args:
        mol (moldesign.Molecule): Molecule to saturate

    Returns:
        moldesign.Molecule: New molecule with all valences saturated
    """
    pbmol = mol_to_pybel(mol)
    pbmol.OBMol.AddHydrogens()
    newmol = pybel_to_mol(pbmol, reorder_atoms_by_residue=True)
    mdt.helpers.assign_unique_hydrogen_names(newmol)
    return newmol


@exports
def mol_to_pybel(mdtmol):
    """ Translate a moldesign molecule object into a pybel molecule object.

    Note:
        The focus is on translating topology and biomolecular structure -
        we don't translate any metadata.

    Args:
        mdtmol (moldesign.Molecule): molecule to translate

    Returns:
        pybel.Molecule: translated molecule
    """
    obmol = ob.OBMol()
    obmol.BeginModify()
    atommap = {}
    resmap = {}
    for atom in mdtmol.atoms:
        obatom = obmol.NewAtom()
        obatom.SetAtomicNum(atom.atnum)
        atommap[atom] = obatom
        pos = atom.position.value_in(u.angstrom)
        obatom.SetVector(*pos)

        if atom.residue and atom.residue not in resmap:
            obres = obmol.NewResidue()
            resmap[atom.residue] = obres
            obres.SetChain(str(atom.chain.pdbname)[0] if atom.chain.pdbname else 'Z')
            obres.SetName(str(atom.residue.pdbname) if atom.residue.pdbname else 'UNL')
            obres.SetNum(str(atom.residue.pdbindex) if atom.residue.pdbindex else 0)
        else:
            obres = resmap[atom.residue]

        obres.AddAtom(obatom)
        obres.SetHetAtom(obatom, not atom.residue.is_standard_residue)
        obres.SetAtomID(obatom, str(atom.name))
        obres.SetSerialNum(obatom,
                           mdt.utils.if_not_none(atom.pdbindex, atom.index+1))

    for atom in mdtmol.bond_graph:
        a1 = atommap[atom]
        for nbr, order in mdtmol.bond_graph[atom].items():
            a2 = atommap[nbr]
            if a1.GetIdx() > a2.GetIdx():
                obmol.AddBond(a1.GetIdx(), a2.GetIdx(), order)
    obmol.EndModify()
    pbmol = pb.Molecule(obmol)

    for atom in atommap:
        idx = atommap[atom].GetIdx()
        obatom = obmol.GetAtom(idx)
        obatom.SetFormalCharge(int(atom.formal_charge.value_in(u.q_e)))
    return pbmol


@exports
def pybel_to_mol(pbmol,
                 reorder_atoms_by_residue=False,
                 primary_structure=True,
                 **kwargs):
    """ Translate a pybel molecule object into a moldesign object.

    Note:
        The focus is on translating topology and biomolecular structure - we don't translate any metadata.

    Args:
        pbmol (pybel.Molecule): molecule to translate
        reorder_atoms_by_residue (bool): change atom order so that all atoms in a residue are stored
            contiguously
        primary_structure (bool): translate primary structure data as well as atomic data
        **kwargs (dict): keyword arguments to  moldesign.Molecule __init__ method

    Returns:
        moldesign.Molecule: translated molecule
    """
    newatom_map = {}
    newresidues = {}
    newchains = {}
    newatoms = moldesign.molecules.atomcollections.AtomList([])
    backup_chain_names = list(string.ascii_uppercase)

    for pybatom in pbmol.atoms:
        obres = pybatom.OBAtom.GetResidue()
        name = obres.GetAtomID(pybatom.OBAtom).strip()

        if pybatom.atomicnum == 67:
            print(("WARNING: openbabel parsed atom serial %d (name:%s) as Holmium; "
                   "correcting to hydrogen. ") % (pybatom.OBAtom.GetIdx(), name))
            atnum = 1

        elif pybatom.atomicnum == 0:
            print("WARNING: openbabel failed to parse atom serial %d (name:%s); guessing %s. " % (
                pybatom.OBAtom.GetIdx(), name, name[0]))
            atnum = moldesign.data.ATOMIC_NUMBERS[name[0]]
        else:
            atnum = pybatom.atomicnum
        mdtatom = mdt.Atom(atnum=atnum, name=name,
                           formal_charge=pybatom.formalcharge * u.q_e,
                           pdbname=name, pdbindex=pybatom.OBAtom.GetIdx())
        newatom_map[pybatom.OBAtom.GetIdx()] = mdtatom
        mdtatom.position = pybatom.coords * u.angstrom

        if primary_structure:
            obres = pybatom.OBAtom.GetResidue()
            resname = obres.GetName()
            residx = obres.GetIdx()
            chain_id = obres.GetChain()
            chain_id_num = obres.GetChainNum()

            if chain_id_num not in newchains:
                # create new chain
                if not mdt.utils.is_printable(chain_id.strip()) or not chain_id.strip():
                    chain_id = backup_chain_names.pop()
                    print('WARNING: assigned name %s to unnamed chain object @ %s' % (
                        chain_id, hex(chain_id_num)))
                chn = mdt.Chain(pdbname=str(chain_id))
                newchains[chain_id_num] = chn
            else:
                chn = newchains[chain_id_num]

            if residx not in newresidues:
                # Create new residue
                pdb_idx = obres.GetNum()
                res = mdt.Residue(pdbname=resname,
                                  pdbindex=pdb_idx)
                newresidues[residx] = res
                chn.add(res)
                res.chain = chn
            else:
                res = newresidues[residx]

            res.add(mdtatom)

        newatoms.append(mdtatom)

    for ibond in range(pbmol.OBMol.NumBonds()):
        obbond = pbmol.OBMol.GetBond(ibond)
        a1 = newatom_map[obbond.GetBeginAtomIdx()]
        a2 = newatom_map[obbond.GetEndAtomIdx()]
        order = obbond.GetBondOrder()
        bond = mdt.Bond(a1, a2)
        bond.order = order

    if reorder_atoms_by_residue and primary_structure:
        resorder = {}
        for atom in newatoms:
            resorder.setdefault(atom.residue, len(resorder))
        newatoms.sort(key=lambda a: resorder[a.residue])

    return mdt.Molecule(newatoms, **kwargs)


def from_smiles(smi, name=None):
    """ Translate a smiles string to a 3D structure.
    This method uses OpenBabel to generate a plausible 3D conformation of the 2D SMILES topology.
    We only use the first result from the conformation generator.

    Args:
        smi (str): smiles string
        name (str): name to assign to molecule (default - the smiles string)

    Returns:
        moldesign.Molecule: the translated molecule
    """
    return _string_to_3d_mol(smi, 'smi', name)


def from_inchi(inchi, name=None):
    """ Translate an INCHI string to a 3D structure.
    This method uses OpenBabel to generate a plausible 3D conformation of the 2D SMILES topology.
    We only use the first result from the conformation generator.

    Args:
        smi (str): smiles string
        name (str): name to assign to molecule (default - the smiles string)

    Returns:
        moldesign.Molecule: the translated molecule
    """
    return _string_to_3d_mol(inchi, 'inchi', name)

@runsremotely(enable=force_remote)
def _string_to_3d_mol(s, fmt, name):
    if name is None: name = s
    pbmol = pb.readstring(fmt, str(s))  # avoid passing unicode by casting to str
    pbmol.addh()
    pbmol.make3D()
    mol = pybel_to_mol(pbmol,
                       name=name,
                       primary_structure=False)
    mdt.helpers.atom_name_check(mol)
    return mol
