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
import cgi
import collections
import functools
import re

import ipywidgets as ipy
import numpy as np

import moldesign as mdt
import moldesign.widgets.logs
from moldesign import compute
from moldesign import units as u
from moldesign.methods import basemethods
from moldesign.utils import if_not_none, DotDict

IMAGE = 'ambertools'

AmberParameters = collections.namedtuple('AmberParameters',
                                         'prmtop inpcrd job')
GAFFParameters = collections.namedtuple('GAFFParameters',
                                        'lib frcmod job')

DNA_FF = {ff: 'leaprc.DNA.%s'%ff for ff in ['OL15', 'bsc1']}
RNA_FF = {ff: 'leaprc.RNA.%s'%ff for ff in ['OL3', 'YIL']}
PROTEIN_FF = {ff: 'leaprc.protein.%s'%ff for ff in
              ['ff14SB',
               'ff14SBonlysc',
               'ff14SB.redq',
               'ff15ipq',
               'ff15ipq-vac',
               'ff03.r1',
               'ff03ua',
               'fb15']}
LIPID_FF = {'lipid14': 'leaprc.lipid14'}
WATER_FF = {ff: 'leaprc.water.%s'%ff for ff in
            ['tip3p', 'tip4pew', 'spce', 'opc']}
CARB_FF = {ff: 'leaprc.%s'%ff for ff in ['GLYCAM_06EPb', 'GLYCAM_06j-1']}
ORGANIC_FF = {ff: 'leaprc.%s'%ff for ff in ['gaff', 'gaff2']}
OTHER_FF = {ff: 'leaprc.%s'%ff for ff in ['phosaa10', 'modrna08', 'xFPchromophores']}
LEAPRCFILES = {}
for ffiles in (DNA_FF, RNA_FF, PROTEIN_FF, LIPID_FF, CARB_FF, WATER_FF, ORGANIC_FF, OTHER_FF):
    LEAPRCFILES.update(ffiles)


class ParameterizationError(Exception): pass


class SQMPotential(basemethods.QMBase):
    DEFAULT_PROPERTIES = ['potential_energy',
                          'electronic_state',
                          'mulliken']
    ALL_PROPERTIES = DEFAULT_PROPERTIES
    THEORIES = 'MNDO MNDO/d AM1 AM1/d PM3 PDDG PDDG/MNDO PDDG/PM3 RM1 PM3CARB1 PM3-MAIS PM6 DFTB'.split()
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

        engine = mdt.compute.default_engine
        image = compute.get_image_path(IMAGE, engine)
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
            if fields[:2] == ['QM','DIPOLE']:
                result['dipole'] = np.array(map(float,fields[2:5])) * u.ureg.debye # TODO CHECK UNITS


def am1_bcc_charges(mol, minsteps=None, wait=True):
    charge = if_not_none(mol.charge, 0)

    engine = mdt.compute.default_engine
    image = compute.get_image_path(IMAGE, engine)
    command = 'antechamber -fi pdb -i mol.pdb -fo mol2 -o out.mol2 -c bcc -an n'
    if charge != 0: command += ' -nc %d' % charge
    if minsteps is not None: command += ' -ek "maxcyc=%d"' % minsteps

    def parse_mol2(job):
        """Callback to complete the job"""
        atom_info = DotDict(job=job)
        lines = iter(job.get_output('out.mol2').read().split('\n'))
        while True:
            line = lines.next()
            fields = line.split()
            if fields[0] == 'ATOM':
                idx = int(fields[1]) - 1
                name = fields[2]
                assert mol.atoms[idx].name == name
                atom_info[mol.atoms[idx]] = DotDict(partialcharge=u.q_e * float(fields[-2]),
                                                    atomtype=fields[-1])
        return atom_info

    job = engine.launch(image,
                          command=command,
                          name="am1-bcc, %s" % mol.name,
                          inputs={'mol.pdb': mol.write(format='pdb')},
                          when_finished=parse_mol2)
    moldesign.widgets.logs.display(job, job.name)
    if not wait: return job()
    else: job.wait()


def build_bdna(sequence,
               image=IMAGE,
               engine=None):
    """
    Uses nucleic acid builder to build a bdna structure.
    This uses very little of NAB's functionality.
    :param sequence: e.g., 'actgaac...'
    :param image: docker image name
    :param engine: compute engine
    :return: moldesign.molecule.Molecule
    """
    from moldesign.structure import converters

    infile = 'molecule m;\nm = bdna( "%s" );\nputpdb( "helix.pdb", m, "-wwpdb");\n' % sequence.lower()
    engine = if_not_none(engine, compute.default_engine)
    imagename = compute.get_image_path(image, engine)

    job = engine.launch(image=imagename,
                        command='nab -o buildbdna build.nab && ./buildbdna',
                        inputs={'build.nab': infile},
                        name='NAB_build_bdna')
    moldesign.widgets.logs.display(job, job.name)
    job.wait()
    mol = converters.read(job.get_output('helix.pdb'), format='pdb')
    mol.name = 'BDNA: %s' % sequence
    return mol


def run_tleap(mol,
              protein='ff14SB',
              dna='OL15',
              rna='OL3',
              carbohydrate='GLYCAM_06j-1',
              lipid='lipid14',
              water='tip3p',
              organic='gaff2',
              off_files=(),
              frcmod_files=(),
              image=IMAGE,
              engine=None):
    """
    Drives tleap to create a prmtop and inpcrd file. Specifically uses the AmberTools 16
    tleap distribution.

    Defaults are as recommended in the ambertools manual.

    Args:
        mol (moldesign.Molecule): Molecule to set up
        protein (str): protein forcefield name (default:ff14SB)
        dna (str): dna forcefield name (default: OL15)
        rna (str): rna forcefield name (default: OL3)
        carbohydrate (str): carbohydrate forcefield name (default: GLYCAM_06j)
        lipid (str): lipid forcefield name (default: lipid14)
        water (str): water forcefield name (default: tip3p)
        organic (str): organic forcefield name (default: gaff2)
        off_files (List[batch.FileContainer]):
        frcmod_files (List[batch.FileContainer]):
        image (str): docker image to use
        engine (batch.EngineBase): compute engine to run tleap on

    References:
        Ambertools Manual, http://ambermd.org/doc12/Amber16.pdf. See page 33 for forcefield
        recommendations.
    """
    # Prepare input for tleap
    leapstr = ['source %s' % LEAPRCFILES[ff] for ff in
               (protein, dna, rna, carbohydrate, lipid, water, organic)]

    for frcmod in frcmod_files:
        fname = frcmod.dumphere()
        leapstr.append('loadamberparam %s' % fname)
    for off in off_files:
        fname = off.dumphere()
        leapstr.append('loadoff %s' % fname)
    leapstr.append('mol = loadpdb input.pdb\n'
                   "check mol\n"
                   "saveamberparm mol output.prmtop output.inpcrd\n"
                   "quit\n")

    # Launch the job
    engine = if_not_none(engine, compute.default_engine)
    imagename = compute.get_image_path(image, engine)
    command = 'tleap -f input.leap'
    inputs = {'input.pdb': mol.write(format='pdb'),
              'input.leap': '\n'.join(leapstr)}

    job = engine.launch(imagename,
                        command,
                        inputs=inputs,
                        name="tleap, %s"% mol.name)
    moldesign.widgets.logs.display(job, "tleap, %s"%mol.name)
    job.wait()

    return job


@functools.wraps(run_tleap)
def assign_forcefield(mol, **kwargs):
    job = run_tleap(mol, **kwargs)
    report = ParamterizationDisplay(job, mol)
    if report.msg:
        moldesign.widgets.logs.display(report, title='ERRORS/WARNINGS')
    if 'output.inpcrd' in job.get_output():
        prmtop = job.get_output('output.prmtop')
        inpcrd = job.get_output('output.inpcrd')
        params = AmberParameters(prmtop, inpcrd, job)
        mol = mdt.read_amber(params.prmtop, params.inpcrd)
        mol.ff.amber_params = params
        return mol
    else:
        raise ParameterizationError('TLeap failed to assign force field parameters for %s' % mol, job)



def get_gaff_parameters(mol, charges, image=IMAGE, engine=None):
    inputs = {}
    cmds = []
    engine = if_not_none(engine, compute.default_engine)
    imagename = compute.get_image_path(image, engine)

    # Add charges to molecule
    inputs['mol.charges'] = '\n'.join(map(str, charges))
    inputs['mol.mol2'] = mol.write(format='mol2')

    # Creates a mol2 file with the passed charges
    cmds.append('antechamber -i mol.mol2 -fi mol2 -o mol_charged.mol2 -fo mol2 -c rc -cf mol.charges')

    # Add missing parameters, write out library and parameter files
    cmds.append('parmchk -i mol_charged.mol2 -f mol2 -o mol.frcmod')

    # Create the lib file
    cmds.append('tleap -f lea.in')
    inputs['leap.in'] = '\n'.join(["source leaprc.%s" % ff,
                                   "source leaprc.gaff",
                                   "LIG = loadmol2 mol_charged.mol2",
                                   "fmod = loadamberparams mol.frcmod",
                                   "check LIG",
                                   "saveoff LIG mol.lib",
                                   "saveamberparm LIG mol.prmtop mol.inpcrd",
                                   "quit\n"])

    # Submit the job and wait
    job = engine.launch(imagename,
                          ' && '.join(cmds),
                          inputs=inputs,
                          name="GAFF assignments" % mol.name)
    moldesign.widgets.logs.display(job, "tleap, %s"%mol.name)
    job.wait()

    param = GAFFParameters(job.get_output('mol.lib'),
                           job.get_output('mol.frcmod'),
                           job)
    return param


class ParamterizationDisplay(ipy.Box):
    ATOMSPEC = re.compile(r'\.R<(\S+) (\d+)>\.A<(\S+) (\d+)>')

    def __init__(self, job, molin, molout=None):
        self.molin = molin
        self.molout = molout
        self.job = job
        self.msg = []
        self.parse_errors(self.job)

        self.listdesc = ipy.HTML('<b>Errors / warnings:</b>')
        self.errorlist = ipy.Select(options=collections.OrderedDict((e.short, e) for e in self.msg))
        self.errmsg = ipy.HTML('-')

        self.viewerpane = self.molin.draw()
        self.viewer = self.viewerpane.children[0]
        self.viewer.ribbon(opacity=0.7)

        if self.errorlist.value is not None:
            self.switch_display({'old': self.errorlist.value, 'new': self.errorlist.value})
        self.errorlist.observe(self.switch_display, 'value')
        children = (ipy.HBox([self.viewerpane, ipy.VBox([self.listdesc, self.errorlist])]),
                    self.errmsg)

        super(ParamterizationDisplay, self).__init__(children=children)

    def switch_display(self, d):
        old = d['old']
        old.unshow(self.viewer)
        self.errmsg.value = '-'
        new = d['new']
        new.show(self.viewer)
        self.errmsg.value = new.desc

    def parse_errors(self, job):
        # TODO: special messages for known problems (e.g. histidine)
        unknown_res = set()
        lineiter = iter(job.stdout.split('\n'))
        reslookup = {str(i + self.molin.residues[0].pdbindex): r for i,r in enumerate(self.molin.residues)}

        def _atom_from_re(s):
            resname, residx, atomname, atomidx = s
            r = reslookup[residx]
            a = r[atomname]
            return a

        def unusual_bond(l):
            atomre1, atomre2 = self.ATOMSPEC.findall(l)
            try:
                a1, a2 = _atom_from_re(atomre1), _atom_from_re(atomre2)
            except KeyError:
                a1 = a2 = None
            r1 = reslookup[atomre1[1]]
            r2 = reslookup[atomre2[1]]
            self.msg.append(UnusualBond(l, (a1, a2), (r1, r2)))

        while True:
            try: line = lineiter.next()
            except StopIteration: break

            fields = line.split()
            if fields[0:2] == ['Unknown','residue:']:
                # EX: "Unknown residue: 3TE   number: 499   type: Terminal/beginning"
                res = self.molin.residues[int(fields[4])]
                self.msg.append(UnknownResidue(line,res))
                unknown_res.add(res)

            elif fields[:4] == 'Warning: Close contact of'.split():
                # EX: "Warning: Close contact of 1.028366 angstroms between .R<DC5 1>.A<HO5' 1> and .R<DC5 81>.A<P 9>"
                unusual_bond(line)

            elif fields[:6] == 'WARNING: There is a bond of'.split():
                # Matches two lines, EX:
                # "WARNING: There is a bond of 34.397700 angstroms between:"
                # "-------  .R<DG 92>.A<O3' 33> and .R<DG 93>.A<P 1>"
                nextline = lineiter.next()
                unusual_bond(line + nextline)

            elif fields[:5] == 'Created a new atom named:'.split():
                # EX: "Created a new atom named: P within residue: .R<DC5 81>"
                residue = reslookup[fields[-1][:-1]]
                if residue in unknown_res: continue  # suppress atoms from an unknown res ...
                atom = residue[fields[5]]
                self.msg.append(UnknownAtom(line, residue, atom))


class ForceFieldMessage(object):
    pass


class UnknownAtom(ForceFieldMessage):
    def __init__(self, message, residue, atom):
        self.message = message
        self.residue = residue
        self.atom = atom
        self.desc = 'ERROR: Atom name "%s" was not found in the "%s" template<br>%s' % (
            self.atom.name, self.residue.resname, self.atom) + '<p>TLeap message:<i>%s</i>' % self.message
        self.short = 'ERR: %s: unknown atom name "%s" for residue %s' % (self.atom,
                                                                         self.atom.name,
                                                                         self.atom.residue.resname)

    def show(self, viewer):
        viewer.licorice(atoms=self.residue.atoms, render=False)
        viewer.vdw(atoms=[self.atom])

    def unshow(self, viewer):
        viewer.ribbon(opacity=0.7)


class MissingAtom(ForceFieldMessage):
    def __init__(self, message, residue, atom):
        self.message = message
        self.residue = residue
        self.atom = atom
        self.desc = 'INFO: Atom %s in %s (chain %s) was added to the system using the "%s" template' % (
            self.atom.name, self.residue.name, self.residue.chain.name, self.residue.resname) + \
                    '<p>TLeap message:<i>%s</i>' % self.message
        self.short = 'INFO: Missing heavy atom %s (index %d)' % (self.atom, self.atom.index)

    def show(self, viewer):
        viewer.licorice(atoms=self.residue.atoms, render=False)
        viewer.vdw(atoms=[self.atom])

    def unshow(self, viewer):
        viewer.ribbon(opacity=0.7)


class UnknownResidue(ForceFieldMessage):
    def __init__(self, message, residue):
        self.message = message
        self.residue = residue
        self._label = None
        self.desc = ('ERROR: Residue type "%s" was not found in residue templates<br>%s' % (self.residue.resname, self.residue)
            + '<p>TLeap message:<i>%s</i>' % self.message)
        self.short = 'ERR: %s: unknown res type "%s"' % (self.residue, self.residue.resname)

    def show(self, viewer):
        viewer.licorice(opacity=1.0, atoms=self.residue.atoms, render=False)
        self._label = viewer.draw_label(position=self.residue.com, text=self.residue.name)

    def unshow(self, viewer):
        viewer.ribbon(opacity=0.7, render=False)
        if self._label: viewer.remove(self._label)
        self._label = None


class UnusualBond(ForceFieldMessage):
    def __init__(self, message, atoms, residues):
        self.message = message
        self.atoms = atoms
        self.residues = residues
        self._has_atoms = (self.atoms[0] is not None) and (self.atoms[1] is not None)
        self._shape = None
        if self._has_atoms:
            self.desc = 'WARNING: Unusual distance between {a1} and {a2}: {d:.3f}'.format(
                d=self.atoms[0].distance(self.atoms[1]), a1=self.atoms[0], a2=self.atoms[1]) \
                        + '<p>TLeap message:<br><i>%s</i>' % cgi.escape(self.message)

            self.short = 'WARN: Unusual dist: {} - {} = ({:.1f})'.format(self.atoms[0],
                                                                   self.atoms[1],
                                                                   self.atoms[0].distance(self.atoms[1]))
        else:
            self.short = 'WARN: Unusual bond - atoms not shown'
            self.desc = 'TLeap message:<br><i>%s</i>' % cgi.escape(self.message)

    def show(self, viewer):
        if self._has_atoms:
            res_opacity = 0.7
        else:
            res_opacity = 0.9

        viewer.licorice(opacity=res_opacity, atoms=self.residues[0], render=False)
        viewer.licorice(opacity=res_opacity, atoms=self.residues[1], render=False)
        if self._has_atoms:
            self._shape = viewer.draw_cylinder(start=self.atoms[0].position,
                                               end=self.atoms[1].position,
                                               radius=0.1 * u.angstrom,
                                               opacity=1.0,
                                               color='red')
        else:
            viewer.render()

    def unshow(self, viewer):
        viewer.ribbon(opacity=0.7, render=False)
        if self._shape: viewer.remove(self._shape)
        self._shape = None

