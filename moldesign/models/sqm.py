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
import numpy as np

import moldesign as mdt
from moldesign import units as u, compute
from moldesign.interfaces.ambertools import IMAGE
from .base import QMBase


class SQMPotential(QMBase):
    DEFAULT_PROPERTIES = ['potential_energy',
                          'wfn',
                          'mulliken']
    ALL_PROPERTIES = DEFAULT_PROPERTIES
    THEORIES = ('MNDO MNDO/d AM1 AM1/d PM3 PDDG PDDG/MNDO PDDG/PM3 RM1 '
               'PM3CARB1 PM3-MAIS PM6 DFTB').split()
    FORCE_UNITS = u.hartree / u.bohr

    def __init__(self, **kwargs):
        super(SQMPotential, self).__init__(**kwargs)

    def calculate(self, requests=None, guess=None, wait=True):
        inputfile = ['SQM input for %s' % self.mol.name,
                     ' &qmmm',
                     " qm_theory='%s' qmcharge=%d, printcharges=1, maxcyc=0" % (
                         self.params.theory, self.get_formal_charge()),
                     '/']
        for atom in self.mol.atoms:
            inputfile.append(' {atom.atnum}  {atom.name}  {pos[0]:.8f} {pos[1]:.8f} {pos[2]:.8f}'.format(
                atom=atom, pos=atom.position.value_in(u.angstrom)))

        engine = mdt.compute.get_engine()
        image = compute.get_image_path(IMAGE)
        job = engine.launch(image,
                              'sqm -i mol.in -o mol.out',
                              inputs={'mol.in': '\n'.join(inputfile)},
                              name="sqm single point, %s" % self.mol.name)
        if not wait: return job

        else: return self._parse_results(job)

    def _parse_results(self, job):
        result = {}
        lines = iter(job.get_output('mol.out').split('\n'))
        while True:
            line = lines.next()
            fields = line.split()
            if fields[:2] == ['QM', 'DIPOLE']:
                # TODO: CHECK UNITS
                result['dipole'] = np.array(map(float,fields[2:5])) * u.ureg.debye