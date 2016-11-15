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
import pyccc

import moldesign as mdt
from moldesign import units as u

from .base import QMBase

IMAGE = 'nwchem'

TASKNAMES = {'rhf': 'scf',
             'uks': 'scf',
             'rks': 'dft',
             'uks': 'dft'}

SCFNAMES = {'rhf': 'rhf',
            'uks': 'uhf',
            'rks': 'rhf',
            'uks': 'uhf'}

class NWChemQM(QMBase):

    def prep(self):
        self._scfname = SCFNAMES[self.params.theory]
        self._taskname = TASKNAMES[self.params.theory]
        self._prepped = True

    def calculate(self, requests=None):
        self.prep()

        grad = (requests is not None and 'forces' in requests)

        job = pyccc.Job(image=mdt.compute.get_image_path(IMAGE),
                        command='nwchem nw.in && get_results.py',
                        inputs=self.make_input_files(grad=grad),
                        when_finished=self.get_results)

        return mdt.compute.run_job(job, _return_result=True)

    def make_input_files(self, grad=False):
        nwin = ['title %s\nstart' % self.mol.name,
                self._header(),
                self._geom_block(self.mol),
                self._basisblock(),
                self._chargeblock(),
                self._theoryblock(),
                self._taskblock(grad)]

        return {'nw.in': '\n'.join(nwin)}

    def get_results(self):
        raise NotImplementedError()

    @staticmethod
    def _geom_block(mol):
        lines = ['geometry units angstrom noautoz noautosym']
        for atom in mol.atoms:
            lines.append('  %s  %24.14e  %24.14e  %24.14e' %
                         (atom.symbol,
                          atom.x.value_in(u.angstrom),
                          atom.y.value_in(u.angstrom),
                          atom.z.value_in(u.angstrom)))
        lines.append('end')
        return '\n'.join(lines)

    def _basisblock(self):
        return 'basis\n  * library %s\nend' % self.params.basis  # TODO: translate names

    def _theoryblock(self):
        lines = [self._taskname,
                 self._multiplicityline(),
                 self._theorylines(),
                 'end'
                 ]
        return '\n'.join(lines)

    def _taskblock(self, grad):
        if grad:
            tasktype = 'gradient'
        else:
            tasktype = 'energy'
        return 'task %s %s' % (self._taskname, tasktype)

    def _multiplicityline(self):
        return 'mult %s' % self.params.multiplicity

    def _chargeblock(self):
        return '\ncharge %s\n' % self.mol.charge.value_in(u.q_e)

    def _header(self):
        return '\nstart mol\n\npermanent_dir ./perm\n'

    def _theorylines(self):
        if self._taskname == 'dft':
            return 'XC %s' % self.params.functional
        else:
            return ''

