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

from .. import data
from ..utils import doc_inherit, exports
from .forcefieldbase import Forcefield

from ..interfaces import ambertools



@exports
class TLeapForcefield(Forcefield):
    """ Amber-type forcefield.

    Assignment routines for these forcefields are backed by ambertools in general (and
    tleap specifically)

    Args:
        fflines (List[str]): commands to load this forcefield
        file_list (Dict[str, pyccc.FileReference]): any files necessary for assigning parameters
             (e.g. lib files, frcmod files)
    """
    TYPE = 'amber'

    def __init__(self, fflines, file_list=None):
        self._fflines = fflines
        self._file_list = file_list if file_list is not None else {}
        super(TLeapForcefield, self).__init__()

    @doc_inherit
    def assign(self, mol, display=True):
        clean_molecule = ambertools._prep_for_tleap(mol)

        job = ambertools.run_tleap(clean_molecule)

        if 'output.inpcrd' in job.get_output():
            prmtop = job.get_output('output.prmtop')
            inpcrd = job.get_output('output.inpcrd')
            params = ambertools.AmberParameters(prmtop, inpcrd, job)
            m = mdt.read_amber(params.prmtop, params.inpcrd)
            newmol = mdt.helpers.restore_topology(m, mol)
            newmol.ff = mdt.forcefields.ForcefieldParams(newmol, params)
        else:
            newmol = None

        errors = ambertools._parse_tleap_errors(job, clean_molecule)

        display.show_parameterization_results(errors, clean_molecule, molout=newmol)

        if newmol is not None:
            return newmol
        else:
            raise ParameterizationError('TLeap failed to assign force field parameters for %s'%mol,
                                        job)

    @doc_inherit
    def add_ff(self, ff):
        self._fflines.extend(ff._fflines)
