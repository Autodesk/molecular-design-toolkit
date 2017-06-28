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
import imp
from ..models import psi4 as psi4_mdt

from moldesign import units as u
from ..utils import exports
from .. import molecules
from .. import units as u


try:
    imp.find_module('psi4')
except (ImportError, OSError) as exc:
    print('Psi4 not installed; using remote docker container')
    force_remote = True
else:
    force_remote = False


@exports
def mdt_to_psi4(mol):
    import psi4
    lines = []
    try:
        assert mol.multiplicity
    except AttributeError:
        pass
    else:
        lines.append('\n')
        lines.append( str(mol.charge) + ' ' + str(mol.multiplicity) )
    
    for atom in mol.atoms:
        x, y, z = atom.position.value_in(u.angstrom)
        lines.append('%s %f %f %f' % (atom.symbol, x, y, z))

    print ('\n'.join(lines))

    geom = psi4.core.Molecule.create_molecule_from_string('\n'.join(lines))
    return geom
#name JSPB June 15
@exports
def mol_name(mol):
    return mol.name

@exports
def Opt_Trajectory(mol, psi4_history):
    import psi4
    my_trajectory=psi4_mdt.capture_history(name_string=mol.name,
                                           molecule=mol,
                                           psi4_history=psi4_history,
                                           model=mol.energy_model)
    return my_trajectory

@exports
def psi4_to_mdt(psi4_molecule):
    import psi4
    
    """

    Args:
        geom (psi4.core.Molecule):

    Returns:

    """
    
    psi4_molecule.update_geometry()
    
    atom_list=[]
    
    for iat in range(psi4_molecule.natom()):
        print(psi4_molecule.symbol(iat))
        atom_list.append(molecules.Atom(psi4_molecule.symbol(iat)))
        atom_list[iat].x = psi4_molecule.x(iat) * u.angstrom
        atom_list[iat].y = psi4_molecule.y(iat) * u.angstrom
        atom_list[iat].z = psi4_molecule.z(iat) * u.angstrom


    my_mol=molecules.molecule.Molecule( atom_list , charge=psi4_molecule.molecular_charge(), multiplicity = psi4_molecule.multiplicity())

    return my_mol




