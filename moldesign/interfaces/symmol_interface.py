from __future__ import print_function, absolute_import, division

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
import re

import numpy as np

import moldesign as mdt
from .. import units as u
from .. import utils
from ..compute import packages


IMAGE = 'symmol'

#@doi('10.1107/S0021889898002180')
def run_symmol(mol, tolerance=0.1 * u.angstrom):
    import fortranformat
    line_writer = fortranformat.FortranRecordWriter('(a6,i2,6f9.5)')

    if mol.num_atoms == 2:
        return _get_twoatom_symmetry(mol)

    infile = ['1.0 1.0 1.0 90.0 90.0 90.0',  # line 1: indicates XYZ coordinates
              # line 2: numbers indicate: mass weighted moment of inertia,
              #         tolerance interpretation, tolerance value,
              #         larger tolerance value (not used)
              '1 0 %f 0.0' % tolerance.value_in(u.angstrom)]
    for atom in mol.atoms:
        infile.append(line_writer.write((atom.element, 1,
                                         atom.x.value_in(u.angstrom),
                                         atom.y.value_in(u.angstrom),
                                         atom.z.value_in(u.angstrom),
                                         0.0, 0.0, 0.0)))
    infile.append('')
    command = 'symmol < sym.in'
    inputs = {'sym.in': '\n'.join(infile)}

    job = packages.symmol.make_job(command=command,
                                   inputs=inputs,
                                   name="symmol, %s" % mol.name)
    job = mdt.compute.run_job(job)

    data = parse_output(mol, job.get_output('symmol.out'))
    _prune_symmetries(data)
    symm = mdt.geom.MolecularSymmetry(
            mol, data.symbol, data.rms,
            orientation=get_aligned_coords(mol, data),
            elems=data.elems,
            _job=job)
    return symm


def _get_twoatom_symmetry(mol):
    """ Symmol doesn't deal with continuous symmetries, so this is hardcoded
    """
    com_pos = mol.positions-mol.com

    ident = mdt.geom.SymmetryElement(mol,
                                     idx=0,
                                     symbol='C1',
                                     matrix=np.identity(3),
                                     csm=0.0*u.angstrom,
                                     max_diff=0.0*u.angstrom)

    # Note: for a continuous symmetry, so the 'matrix' should really be an infinitesimal transform.
    # This doesn't come up a lot practically, so we just put a small generator here instead.
    # The rotation is through an irrational angle so that it does in fact generate the full
    # symmetry group
    transmat = mdt.external.transformations.rotation_matrix(1/np.sqrt(20), com_pos[0])
    axis_rot = mdt.geom.SymmetryElement(mol,
                                        idx=1,
                                        symbol='Cinf_v',
                                        matrix=transmat[:3,:3],
                                        csm=0.0*u.angstrom,
                                        max_diff=0.0*u.angstrom)

    elems = [ident, axis_rot]

    # This is for the mirror plane / inversion center between the two atoms
    if mol.atoms[0].atnum == mol.atoms[1].atnum and mol.atoms[0].mass == mol.atoms[1].mass:
        reflmat = mdt.external.transformations.reflection_matrix([0,0,0], com_pos[0])
        elems.append(mdt.geom.SymmetryElement(mol,
                                              idx=2,
                                              symbol='Cs',
                                              matrix=reflmat[:3,:3],
                                              csm=0.0*u.default.length,
                                              max_diff=0.0*u.default.length))

        elems.append(mdt.geom.SymmetryElement(mol,
                                              idx=2,
                                              symbol='Ci',
                                              matrix=-1 * np.identity(3),
                                              csm=0.0*u.default.length,
                                              max_diff=0.0*u.default.length))

        term_symbol = 'Dinf_h'
    else:
        term_symbol = axis_rot.symbol

    symm = mdt.geom.MolecularSymmetry(mol, term_symbol, 0.0 * u.default.length,
                                      orientation=com_pos,
                                      elems=elems)
    return symm


def _prune_symmetries(data):
    """ Remove identical symmetries
    """
    found = {}
    for elem in data.elems:
        found.setdefault(elem.symbol, [])
        for otherelem in found[elem.symbol]:
            if (np.abs(elem.matrix.T - otherelem.matrix) < 1e-11).all():
                break
        else:
            found[elem.symbol].append(elem)

    newelems = []
    for val in found.values():
        newelems.extend(val)

    data.elems = newelems


MATRIXSTRING = 'ORTHOGONALIZATION MATRIX'.split()
ELEMENTSTRING = 'Symmetry element its CSM and Max.Diff.  Symmetry element its CSM and Max.Diff.'.split()
TRANSFORMATIONSTRING = 'SYMMETRY GROUP MATRICES'.split()
NOSYMM = 'NO SYMMETRY EXISTS WITHIN GIVEN TOLERANCE'.split()
ELEMPARSER = re.compile('(\d+)\) \[(...)\]\s+(\S+)\s+([\-0-9\.]+)\s+([\-0-9\.]+)')
# this regex parses '1) [E  ]  x,y,z       0.0000  0.0000' -> [1, 'E  ', 'x,y,z','0.0000','0.0000']
MATRIXPARSER = re.compile('(\d+)\s+CSM =\s+([\d\.]+)\s+MAX. DIFF. \(Angstrom\)=([\d\.]+)\s+TYPE (\S+)')
# this parses '  4 CSM =   0.06     MAX. DIFF. (Angstrom)=0.0545     TYPE C3' -> [4, 0.06, 0.545, C3]
def parse_output(mol, outfile):
    lines = iter(outfile)
    data = utils.DotDict()
    while True:
        l = next(lines)
        fields = l.split()
        if fields == NOSYMM:
            data.symbol = 'C1'
            data.rms = data.cms = 0.0 * u.angstrom
            data.elems = []
            data.orthmat = np.identity(3)
            return data

        elif fields == MATRIXSTRING:  # get coordinates along principal axes
            data.orthmat = np.zeros((3, 3))
            for i in range(3):
                data.orthmat[i] = list(map(float, next(lines).split()))

        elif fields[:2] == 'Schoenflies symbol'.split():
            data.symbol = fields[3]
            data.csm = float(fields[6]) * u.angstrom
            data.rms = float(fields[-1]) * u.angstrom

        elif fields == ELEMENTSTRING:
            data.elems = []
            while True:
                try:
                    l = next(lines)
                except StopIteration:
                    break
                if l.strip() == '': break
                parsed = ELEMPARSER.findall(l)
                for p in parsed:
                    elem = mdt.geom.SymmetryElement(mol,
                                                    idx=int(p[0])-1,
                                                    symbol=p[1].strip(),
                                                    matrix=_string_to_matrix(p[2]),
                                                    csm=float(p[3]),
                                                    max_diff=float(p[4]) * u.angstrom)
                    if elem.symbol == 'E': elem.symbol = 'C1'
                    data.elems.append(elem)
            break

        elif fields == TRANSFORMATIONSTRING:
            data.elems = []
            l = next(lines)

            while True:
                while l.strip() == '':
                    try: l = next(lines)
                    except StopIteration: return data
                eleminfo = MATRIXPARSER.findall(l)

                if not eleminfo:
                    return data  # we're done
                assert len(eleminfo) == 1
                info = eleminfo[0]

                matrix = np.zeros((3, 3))
                for i in range(3):
                    l = next(lines)
                    matrix[i, :] = list(map(float, l.split()))

                e = mdt.geom.SymmetryElement(mol,
                                             matrix=matrix,
                                             idx=int(info[0])-1,
                                             csm=float(info[1]) * u.angstrom,
                                             max_diff=float(info[2]) * u.angstrom,
                                             symbol=info[3])
                if e.symbol == 'E':
                    e.symbol = 'C1'

                e.matrix = matrix
                data.elems.append(e)
                l = next(lines)

    return data


DIMNUMS = {'x': 0, 'y': 1, 'z': 2}
TRANSFORM_PARSER = re.compile('([+\-]?)([0-9\.]*)([xyz])')
# this regex transforms '3z-5.3x+y' -> [('','3','z'),('-','5.3','x'),('+','','y')]
def _string_to_matrix(string):
    """
    Symmol often returns a symmetry operation as something like "+x,-z,+y"
    This means that the x axis is mapped onto x, y axis is mapped onto -z, and +y is mapped onto z.
    We translate this into a transformation matrix here.
    :param string: A string representing axis mapping, of the form 'x,-z,x+y'
    :return: 3x3 transformation matrix
    """
    mat = []
    for dim in string.split(','):
        row = np.zeros(3)
        components = TRANSFORM_PARSER.findall(dim)
        for sign, factor, dimname in components:
            if factor.strip() == '': factor = '1'
            idim = DIMNUMS[dimname]
            row[idim] = float(sign + factor)
        row = row / np.sqrt(row.dot(row))  # normalize
        mat.append(row)
    return np.array(mat)


def get_aligned_coords(mol, data):
    com = mol.com
    centerpos = mol.positions - com
    orthcoords = (centerpos.T.ldot(data.orthmat)).T
    return orthcoords


