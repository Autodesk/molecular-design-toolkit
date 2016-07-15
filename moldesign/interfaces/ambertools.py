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
import collections

import moldesign as mdt
import pyccc
from moldesign import compute, uibase, utils
from moldesign import units as u

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


def am1_bcc_charges(mol, minsteps=None, wait=True):
    """ Doesn't work yet ..."""
    charge = utils.if_not_none(mol.charge, 0)

    engine = mdt.compute.default_engine
    image = compute.get_image_path(IMAGE, engine)
    command = 'antechamber -fi pdb -i mol.pdb -fo mol2 -o out.mol2 -c bcc -an n'
    if charge != 0: command += ' -nc %d' % charge
    if minsteps is not None: command += ' -ek "maxcyc=%d"' % minsteps

    def parse_mol2(job):
        """Callback to complete the job"""
        atom_info = utils.DotDict(job=job)
        lines = iter(job.get_output('out.mol2').read().split('\n'))
        while True:
            line = lines.next()
            fields = line.split()
            if fields[0] == 'ATOM':
                idx = int(fields[1]) - 1
                name = fields[2]
                assert mol.atoms[idx].name == name
                atom_info[mol.atoms[idx]] = utils.DotDict(partialcharge=u.q_e * float(fields[-2]),
                                                    atomtype=fields[-1])
        return atom_info

    job = engine.launch(image,
                        command=command,
                        name="am1-bcc, %s" % mol.name,
                        inputs={'mol.pdb': mol.write(format='pdb')},
                        when_finished=parse_mol2)
    uibase.display_log(job.get_display_object(), job.name)
    if not wait: return job()
    else: job.wait()


@utils.kwargs_from(mdt.compute.run_job)
def build_bdna(sequence, **kwargs):
    """ Uses Ambertools' Nucleic Acid Builder to build a 3D double-helix B-DNA structure.

    Args:
        sequence (str): DNA sequence for one of the strands (a complementary sequence will
            automatically be created)
        **kwargs: arguments for :meth:`compute.run_job`

    Returns:
        moldesign.Molecule: B-DNA double helix
    """
    infile = 'molecule m;\nm = bdna( "%s" );\nputpdb( "helix.pdb", m, "-wwpdb");\n'% \
             sequence.lower()

    def finish_job(job):
        mol = mdt.read(job.get_output('helix.pdb'), format='pdb')
        mol.name = 'BDNA: %s' % sequence
        return mol

    job = pyccc.Job(command='nab -o buildbdna build.nab && ./buildbdna',
                    image=mdt.compute.get_image_path(IMAGE),
                    inputs={'build.nab': infile},
                    name='NAB_build_bdna',
                    when_finished=finish_job)

    return mdt.compute.run_job(job, _return_result=True, **kwargs)


# TODO: use a single force field specification object rather than 20 kwargs
@utils.kwargs_from(compute.run_job)
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
              **kwargs):
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
        **kwargs: keyword arguments to :meth:`compute.run_job`

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
    inputs = {'input.pdb': mol.write(format='pdb'),
              'input.leap': '\n'.join(leapstr)}

    job = pyccc.Job(image=compute.get_image_path(IMAGE),
                    command='tleap -f input.leap',
                    inputs=inputs,
                    name="tleap, %s" % mol.name)

    return compute.run_job(job, **kwargs)


@mdt.utils.args_from(run_tleap)
def assign_forcefield(mol, **kwargs):
    """ see run_tleap docstring """
    from moldesign.widgets.parameterization import ParameterizationDisplay
    job = run_tleap(mol, **kwargs)

    if 'output.inpcrd' in job.get_output():
        prmtop = job.get_output('output.prmtop')
        inpcrd = job.get_output('output.inpcrd')
        params = AmberParameters(prmtop, inpcrd, job)
        newmol = mdt.read_amber(params.prmtop, params.inpcrd)
        newmol.ff.amber_params = params
    else:
        newmol = None

    report = ParameterizationDisplay(job, mol, molout=newmol)
    uibase.display_log(report, title='ERRORS/WARNINGS', show=True)

    if newmol is not None:
        return newmol
    else:
        raise ParameterizationError('TLeap failed to assign force field parameters for %s' % mol, job)


def get_gaff_parameters(mol, charges, image=IMAGE, engine=None):
    """ Doesn't work yet"""
    inputs = {}

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
    uibase.display_log(job.get_display_object(), "tleap, %s"%mol.name)
    job.wait()

    param = GAFFParameters(job.get_output('mol.lib'),
                           job.get_output('mol.frcmod'),
                           job)
    return param


