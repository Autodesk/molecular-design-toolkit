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
import itertools
import json
import numpy as np

import pyccc

import moldesign as mdt

from moldesign.utils import exports
from moldesign import units as u

from .base import QMBase, QMMMBase
from .jsonmodel import JsonModelBase


@exports
class NWChemQM(JsonModelBase, QMBase):
    """ Interface with NWChem package (QM only)

    Note:
        This is the first interface based on our new wrapping strategy. This is slightly hacked,
        but has the potential to become very general; very few things here are NWChem-specific
    """
    IMAGE = 'nwchem'
    MODELNAME = 'nwchem'
    DEFAULT_PROPERTIES = ['potential_energy']
    ALL_PROPERTIES = DEFAULT_PROPERTIES + 'forces dipole esp'.split()

    def _get_inputfiles(self):
        return {'input.xyz': self.mol.write(format='xyz')}

    def write_constraints(self, parameters):
        parameters['constraints'] = []
        for constraint in self.mol.constraints:

            if not constraint.satisfied():  # TODO: factor this out into NWChem subclass
                raise ValueError('Constraints must be satisfied before passing to NWChem %s'
                                 % constraint)

            cjson = {'type': constraint['desc'],
                     'value': constraint['value'].to_json()}
            if cjson['type'] == 'position':
                cjson['atomIdx'] = constraint.atom.index
            elif cjson['type'] == 'distance':
                cjson['atomIdx1'] = constraint.a1.index
                cjson['atomIdx2'] = constraint.a2.index


@exports
class NWChemQMMM(NWChemQM):
    """ Interface with NWChem package for QM/MM only. Note that this is currently only set up for
    optimizations, and only of the QM geometry - the MM region cannot move.
    """
    IMAGE = 'nwchem'
    MODELNAME = 'nwchem_qmmmm'
    DEFAULT_PROPERTIES = ['potential_energy', 'forces', 'esp']
    ALL_PROPERTIES = DEFAULT_PROPERTIES
    RUNNER = 'runqmmm.py'
    PARSER = 'getresults.py'
    PARAMETERS = NWChemQM.PARAMETERS + [mdt.parameters.Parameter('qm_atom_indices')]

    def _get_inputfiles(self):
        crdparamfile = self._makecrdparamfile()
        return {'nwchem.crdparams': crdparamfile}

    def _makecrdparamfile(self, include_mmterms=True):
        pmdparms = self.mol.ff.parmed_obj

        lines = ['qm']
        qmatom_idxes = set(self.params.qm_atom_indices)

        def crosses_boundary(*atoms):
            total = 0
            inqm = 0
            for atom in atoms:
                total += 1
                if atom.idx in qmatom_idxes:
                    inqm += 1
            return inqm != 0 and inqm < total

        # write QM atoms
        for atomidx in sorted(self.params.qm_atom_indices):
            atom = self.mol.atoms[atomidx]
            x, y, z = atom.position.value_in(u.angstrom)
            lines.append('  %d %s   %24.14f  %24.14f  %24.14f' %
                         (atom.index+1, atom.element, x, y, z))
            qmatom_idxes.add(atom.index)

        # write MM atoms
        lines.append('end\n\nmm')
        for atom, patm in zip(self.mol.atoms, pmdparms.atoms):
            assert atom.index == patm.idx
            assert atom.atnum == patm.element
            x, y, z = atom.position.value_in(u.angstrom)
            if atom.index not in qmatom_idxes:
                lines.append('  %d %s   %24.14f  %24.14f  %24.14f  %24.14f'
                             %(atom.index+1, atom.element, x, y, z, patm.charge))

        lines.append('end\n\nbond\n#      i       j       k_ij             r0')
        if include_mmterms:
            for term in pmdparms.bonds:
                if crosses_boundary(term.atom1, term.atom2):
                    lines.append('  %d    %d   %20.10f    %20.10f' %
                                 (term.atom1.idx+1, term.atom2.idx+1, term.type.k, term.type.req))

        lines.append('end\n\nangle\n#      i       j       k      k_ijk           theta0')
        if include_mmterms:
            for term in pmdparms.angles:
                if crosses_boundary(term.atom1, term.atom2, term.atom3):
                    lines.append('  %d    %d   %d   %20.10f    %20.10f' %
                                 (term.atom1.idx+1, term.atom2.idx+1, term.atom3.idx+1,
                                  term.type.k, term.type.theteq))

        lines.append('end\n\ndihedral\n'
                     '#      i       j       k       l      k_ijkl       periodicity        phase')
        if include_mmterms:
            for term in pmdparms.dihedrals:
                if crosses_boundary(term.atom1, term.atom2, term.atom3, term.atom4):
                    lines.append('  %d    %d   %d   %d   %20.10f   %d    %20.10f' %
                                 (term.atom1.idx+1, term.atom2.idx+1, term.atom3.idx+1, term.atom4.idx+1,
                                  term.type.phi_k, term.type.per, term.type.phase))

        lines.append('end\n\nvdw\n'
                     '#      i       j    A_coeff         B_coeff')
        if include_mmterms:
            for atom1idx, atom2 in itertools.product(qmatom_idxes, pmdparms.atoms):
                if atom2.idx in qmatom_idxes: continue
                atom1 = pmdparms.atoms[atom1idx]
                epsilon = np.sqrt(atom1.epsilon * atom2.epsilon)
                if epsilon == 0: continue
                sigma = (atom1.sigma + atom2.sigma) / 2.0
                lj_a = 4.0 * epsilon * sigma**12
                lj_b = 4.0 * epsilon * sigma**6
                lines.append('  %d    %d  %20.10e   %20.10e' %
                             (atom1.idx+1, atom2.idx+1, lj_a, lj_b))

        lines.append('end\n\nscaled_vdw\n'
                     '#      i       j    A_coeff         B_coeff          one_scnb')
        # TODO: this part
        lines.append('end')

        return pyccc.files.StringContainer('\n'.join(lines))


    def finish_min(self, job):
        traj = mdt.Trajectory(self.mol)
        traj.new_frame()

        results = json.loads(job.get_output('results.json').read())
        new_state = self._json_to_quantities(results['states'][0])

        for iatom in sorted(self.params.qm_indices):
            for position in new_state['positions']:
                self.mol.atoms[iatom].position = position

        properties = self._process_results(results)
        properties.positions = self.mol.positions
        self.mol.properties = properties
        traj.new_frame()
        return traj

