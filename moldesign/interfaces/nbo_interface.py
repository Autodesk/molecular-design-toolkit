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
import moldesign.molecules.bonds
from moldesign import units as u
from moldesign import utils
from moldesign.utils import DotDict


SIGMA_UTF = u"\u03C3"
PI_UTF = u"\u03C0"


def run_nbo(mol, requests=('nlmo', 'nbo'),
            image='nbo',
            engine=None):
    wfn = mol.wfn
    inputs = {'in.47': make_nbo_input_file(mol, requests)}
    command = 'gennbo.i4.exe in.47'
    engine = utils.if_not_none(engine, mdt.compute.config.get_engine())
    imagename = mdt.compute.get_image_path(image)
    job = engine.launch(imagename,
                          command,
                          inputs=inputs,
                          name="nbo, %s" % mol.name)
    moldesign.uibase.display_log(job.get_display_object(), "nbo, %s"%mol.name)

    job.wait()
    parsed_data = parse_nbo(job.get_output('FILE.10'),
                            len(mol.wfn.aobasis))

    for orbtype, data in parsed_data.iteritems():
        if orbtype[0] == 'P':  # copy data from the orthogonal orbitals
            orthdata = parsed_data[orbtype[1:]]
            for key in 'bond_names iatom jatom stars bondnums num_bonded_atoms'.split():
                data[key] = orthdata[key]
            data.occupations = [None for orb in data.coeffs]

        add_orbitals(mol, wfn, data, orbtype)
    wfn._nbo_job = job


def add_orbitals(mol, wfn, orbdata, orbtype):
    orbs = []
    for i in xrange(len(orbdata.coeffs)):
        bond = None
        atoms = [mol.atoms[orbdata.iatom[i] - 1]]
        if orbdata.bond_names[i] == 'RY':
            bname = '%s Ryd*' % atoms[0].name
            nbotype = 'rydberg'
            utf_name = bname
        elif orbdata.bond_names[i] == 'LP':
            bname = '%s lone pair' % atoms[0].name
            nbotype = 'lone pair'
            utf_name = bname
        elif orbdata.bond_names[i] == 'LV':
            bname = '%s lone vacancy' % atoms[0].name
            nbotype = 'lone vacancy'
            utf_name = bname
        elif orbdata.num_bonded_atoms[i] == 1:
            bname = '%s Core' % atoms[0].name
            nbotype = 'core'
            utf_name = bname
        else:
            atoms.append(mol.atoms[orbdata.jatom[i] - 1])
            bond = moldesign.molecules.bonds.Bond(*atoms)
            if orbdata.bondnums[i] == 1:  # THIS IS NOT CORRECT
                nbotype = 'sigma'
                utf_type = SIGMA_UTF
            else:
                nbotype = 'pi'
                utf_type = PI_UTF
            bname = '%s%s (%s - %s)' % (nbotype, orbdata.stars[i],
                                        atoms[0].name, atoms[1].name)
            utf_name = '%s%s (%s - %s)' % (utf_type, orbdata.stars[i],
                                           atoms[0].name, atoms[1].name)
        name = '%s %s' % (bname, orbtype)

        orbs.append(moldesign.methods.orbitals.Orbital(orbdata.coeffs[i],
                                                       wfn=wfn, occupation=orbdata.occupations[i],
                                                       atoms=atoms, name=name,
                                                       nbotype=nbotype,
                                                       bond=bond,
                                                       unicode_name=utf_name,
                                                       _data=orbdata))
    return wfn.add_orbitals(orbs, orbtype=orbtype)


def make_nbo_input_file(mol, requests):
    """
    :param mol:
    :type mol: moldesign.molecules.Molecule
    :return:
    """
    # Units: angstroms, hartrees
    wfn = mol.wfn
    orbs = wfn.molecular_orbitals
    nbofile = []
    # TODO: check for open shell wfn (OPEN keyword)
    # TODO: check normalization, orthogonalization
    nbofile.append(" $GENNBO BODM NATOMS=%d NBAS=%d $END" %
                   (mol.num_atoms, len(wfn.aobasis)))

    commands = ['NBOSUM']
    for r in requests:
        commands.append('AO%s=W10' % r.upper())
        if r[0] != 'P': commands.append('%s' % r.upper())
    nbofile.append('$NBO %s $END' % (' '.join(commands)))

    nbofile.append("$COORD\n %s" % mol.name)
    for iatom, atom in enumerate(mol.atoms):
        #TODO: deal with pseudopotential electrons
        x, y, z = map(repr, atom.position.value_in(u.angstrom))
        nbofile.append("%d %d %s %s %s" % (atom.atnum, atom.atnum,
                                                    x, y, z))
    nbofile.append("$END")

    nbofile.append("$BASIS")
    nbofile.append(' CENTER = ' +
                   ' '.join(str(1+bfn.atom.index) for bfn in wfn.aobasis))
    nbofile.append(" LABEL = " +
                   ' '.join(str(AOLABELS[bfn.orbtype]) for bfn in wfn.aobasis))
    nbofile.append('$END')

    #TODO: deal with CI wavefunctions ($WF keyword)
    nbofile.append('$OVERLAP')
    append_matrix(nbofile, wfn.aobasis.overlaps)
    nbofile.append('$END')

    nbofile.append('$DENSITY')
    append_matrix(nbofile, wfn.density_matrix)
    nbofile.append('$END')
    return '\n '.join(nbofile)


def parse_nbo(f, nbasis):
    lines = f.__iter__()
    parsed = {}
    while True:
        try:
            l = lines.next()
        except StopIteration:
            break

        fields = l.split()
        if fields[1:5] == 'in the AO basis:'.split():
            orbname = fields[0]
            assert orbname[-1] == 's'
            orbname = orbname[:-1]
            lines.next()
            if orbname[0] == 'P':  # these are pre-orthogonal orbitals, it only prints the coefficients
                coeffs = _parse_wrapped_matrix(lines, nbasis)
                parsed[orbname] = DotDict(coeffs=np.array(coeffs))
            else:  # there's more complete information available
                parsed[orbname] = read_orbital_set(lines, nbasis)
    return parsed


def read_orbital_set(lineiter, nbasis):
    # First, get the actual matrix
    mat = _parse_wrapped_matrix(lineiter, nbasis)

    # First, occupation numbers
    occupations = map(float,_get_wrapped_separated_vals(lineiter, nbasis))

    # Next, a line of things that always appear to be ones (for spin orbitals maybe?)
    oneline = _get_wrapped_separated_vals(lineiter, nbasis)
    for x in oneline: assert x == '1'

    # next is number of atoms involved in the bond
    num_bonded_atoms = map(int, _get_wrapped_separated_vals(lineiter, nbasis))
    bond_names = _get_wrapped_separated_vals(lineiter, nbasis)

    # Next indicates whether real or virtual
    stars = _get_wrapped_column_vals(lineiter, nbasis)
    for s in stars: assert (s == '' or s == '*')

    # number of bonds between this pair of atoms
    bondnums = map(int, _get_wrapped_separated_vals(lineiter, nbasis))

    # first atom index (1-based)
    iatom = map(int, _get_wrapped_separated_vals(lineiter, nbasis))
    jatom = map(int, _get_wrapped_separated_vals(lineiter, nbasis))

    # The rest appears to be 0 most of the time ...

    return DotDict(coeffs=np.array(mat),
                   iatom=iatom, jatom=jatom, bondnums=bondnums,
                   bond_names=bond_names,
                   num_bonded_atoms=num_bonded_atoms,
                   stars=stars, occupations=occupations)


def _parse_wrapped_matrix(lineiter, nbasis):
    mat = []
    for i in xrange(nbasis):
        currline = map(float, _get_wrapped_separated_vals(lineiter, nbasis))
        assert len(currline) == nbasis
        mat.append(currline)
    return mat


def _get_wrapped_separated_vals(lineiter, nbasis):
    vals = []
    while True:
        l = lineiter.next()
        vals.extend(l.split())
        if len(vals) == nbasis:
            break
        assert len(vals) < nbasis
    return vals


def _get_wrapped_column_vals(lineiter, nbasis):
    vals = []
    while True:
        l = lineiter.next()[1:]
        lenl = len(l)
        for i in xrange(20):
            if lenl <= 3*i + 1: break
            vals.append(l[3*i: 3*i + 3].strip())
        if len(vals) == nbasis:
            break
        assert len(vals) < nbasis
    return vals


def append_matrix(l, mat):
    for row in mat:
        icol = 0
        while icol < len(row):
            l.append('   ' + ' '.join(map(repr, row[icol:icol + 6])))
            icol += 6


AOLABELS = {'s': 1, 'px': 101, 'py': 102, 'pz': 103,
            "dxx": 201, "dxy": 202, "dxz": 203, "dyy": 204, "dyz": 205, "dzz": 206,
            "fxxx": 301, "fxxy": 302, "fxxz": 303, "fxyy": 304, "fxyz": 305,
            "fxzz": 306, "fyyy": 307, "fyyz": 308, "fyzz": 309, "fzzz": 310,
            "gxxxx": 401, "gxxxy": 402, "gxxxz": 403, "gxxyy": 404, "gxxyz": 405,
            "gxxzz": 406, "gxyyy": 407, "gxyyz": 408, "gxyzz": 409, "gxzzz": 410,
            "gyyyy": 411, "gyyyz": 412, "gyyzz": 413, "gyzzz": 414, "gzzzz": 415,  # end of cartesian
            # start of spherical:
            'p(x)': 151, 'p(y)': 152, 'p(z)': 153,
            "d(xy)": 251, "d(xz)": 252, "d(yz)": 253, "d(x2-y2)": 254, "d(z2)": 255,
            "f(z(5z2-3r2))": 351, "f(x(5z2-r2))": 352, "f(y(5z2-r2))": 353, "f(z(x2-y2))": 354, "f(xyz)": 355,
            "f(x(x2-3y2))": 356, "f(y(3x2-y2))": 357}
