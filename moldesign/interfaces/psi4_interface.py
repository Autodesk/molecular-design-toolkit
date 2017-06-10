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

from moldesign import units as u
from ..utils import exports


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
    for atom in mol.atoms:
        x, y, z = atom.position.value_in(u.angstrom)
        lines.append('%s %f %f %f' % (atom.symbol, x, y, z))

    geom = psi4.geometry('\n'.join(lines))
    return geom

@exports
def psi4_to_mdt(geom):
    """

    Args:
        geom (psi4.core.Molecule):

    Returns:

    """
    geom.update_coords()
    coords = geom.geometry().to_array() * u.bohr




