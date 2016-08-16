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
ExtraAmberParameters = collections.namedtuple('GAFFParameters',
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


class ParameterizationError(Exception):
    pass


@utils.kwargs_from(mdt.compute.run_job)
def calc_am1_bcc_charges(mol, **kwargs):
    """Calculate am1 bcc charges

    Args:
        mol (moldesign.Molecule): assign partial charges to this molecule
            (they will be stored at ``mol.properties['am1-bcc']``)

    Note:
        This will implicity run an AM1 energy minimization before calculating the final
        partial charges. For more control over this process, use the
        ``moldesign.models.SQMPotential`` energy model to calculate the charges.

    Returns:
        Mapping[moldesign.Atom, units.Scalar[charge]]: AM1-BCC partial charges on each atom
    """
    return _antechamber_calc_charges(mol, 'bcc', 'am1-bcc', kwargs)


@utils.kwargs_from(mdt.compute.run_job)
def calc_gasteiger_charges(mol, **kwargs):
    """Calculate gasteiger charges

    Args:
        mol (moldesign.Molecule): assign partial charges to this molecule

    Returns:
        Mapping[moldesign.Atom, units.Scalar[charge]]: gasteiger partial charges on each atom
             (they will be stored at ``mol.properties['gasteiger']``)
    """
    return _antechamber_calc_charges(mol, 'gas', 'gasteiger', kwargs)


def _antechamber_calc_charges(mol, ambname, chargename, kwargs):
    charge = utils.if_not_none(mol.charge, 0)
    command = 'antechamber -fi pdb -i mol.pdb -fo mol2 -o out.mol2 -c %s -an n'%ambname
    if charge != 0:
        command += ' -nc %d'%charge.value_in(u.q_e)

    def finish_job(job):
        """Callback to complete the job"""
        lines = iter(job.get_output('out.mol2').read().split('\n'))
        charges = utils.DotDict(type='atomic')

        line = lines.next()
        while line.strip()[:len('@<TRIPOS>ATOM')] != '@<TRIPOS>ATOM':
            line = lines.next()

        line = lines.next()
        while line.strip()[:len('@<TRIPOS>BOND')] != '@<TRIPOS>BOND':
            fields = line.split()
            idx = int(fields[0])-1
            assert mol.atoms[idx].name == fields[1]
            charges[mol.atoms[idx]] = u.q_e*float(fields[-1])
            line = lines.next()

        mol.properties[chargename] = charges
        return charges

    job = pyccc.Job(image=mdt.compute.get_image_path(IMAGE),
                    command=command,
                    name="%s, %s"%(chargename, mol.name),
                    inputs={'mol.pdb': mol.write(format='pdb')},
                    when_finished=finish_job)
    return compute.run_job(job, _return_result=True, **kwargs)


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


@utils.kwargs_from(compute.run_job)
def run_tleap(mol, forcefields=None, parameters=None, **kwargs):
    """
    Drives tleap to create a prmtop and inpcrd file. Specifically uses the AmberTools 16
    tleap distribution.

    Defaults are as recommended in the ambertools manual.

    Args:
        mol (moldesign.Molecule): Molecule to set up
        forcefields (List[str]): list of the names of forcefields to use
            (see AmberTools manual for descriptions)
        parameters (List[ExtraAmberParameters]): (optional) list of amber parameters
            for non-standard residues
        **kwargs: keyword arguments to :meth:`compute.run_job`

    References:
        Ambertools Manual, http://ambermd.org/doc12/Amber16.pdf. See page 33 for forcefield
        recommendations.
    """
    # Prepare input for tleap
    if forcefields is None: forcefields = mdt.forcefields.ffdefaults.values()
    leapstr = ['source %s' % LEAPRCFILES[ff] for ff in forcefields]
    inputs = {'input.pdb': mdt.helpers.insert_ter_records(mol, mol.write(format='pdb'))}

    if parameters:
        if isinstance(parameters, ExtraAmberParameters):
            parameters = [parameters]
        for ipmtr, p in enumerate(parameters):
            frcname = 'res%d.frcmod' % ipmtr
            libname = 'res%d.lib' % ipmtr
            inputs[frcname] = p.frcmod
            inputs[libname] = p.lib
            leapstr.append('loadamberparam %s' % frcname)
            leapstr.append('loadoff %s' % libname)

    leapstr.append('mol = loadpdb input.pdb\n'
                   "check mol\n"
                   "saveamberparm mol output.prmtop output.inpcrd\n"
                   "savepdb mol output.pdb\n"
                   "quit\n")

    inputs['input.leap'] = '\n'.join(leapstr)

    job = pyccc.Job(image=compute.get_image_path(IMAGE),
                    command='tleap -f input.leap',
                    inputs=inputs,
                    name="tleap, %s" % mol.name)

    return compute.run_job(job, **kwargs)


@utils.args_from(run_tleap)
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


@utils.kwargs_from(mdt.compute.run_job)
def parameterize(mol, charges='esp', ffname='gaff2', **kwargs):
    """Parameterize ``mol``, typically using GAFF parameters.

    Note: this requires partial charges to alread

    Args:
        mol (moldesign.Molecule):
        charges (str): what kind of partial charges to use? These will be taken from
           ``mol.properties.[charges name]``; typical values will be 'esp', 'mulliken',
           'am1-bcc', etc. Use 'zero' to set all charges to 0 (for QM/MM and testing)
           'am1-bcc' will be calculated automatically if not available.
        ffname (str): Name of the gaff-like forcefield file (default: gaff2)

    Returns:
        ExtraAmberParameters: Parameters for the molecule
    """
    assert mol.num_residues == 1
    resname = mol.residues[0].resname

    if charges == 'am1-bcc' and 'am1-bcc' not in mol.properties:
        calc_am1_bcc_charges(mol)
    elif charges == 'gasteiger' and 'gasteiger' not in mol.properties:
        calc_gasteiger_charges(mol)

    if charges == 'zero':
        charges = [0.0 for atom in mol.atoms]
    else:
        charges = u.array([mol.properties[charges][atom] for atom in mol.atoms])
        if not charges.dimensionless:  # implicitly convert floats to fundamental charge units
            charges = charges.to(u.q_e).magnitude

    inputs = {'mol.mol2': mol.write(format='mol2'),
              'mol.charges': '\n'.join(map(str, charges))}

    cmds = ['antechamber -i mol.mol2 -fi mol2 -o mol_charged.mol2 -fo mol2 -c rc -cf mol.charges',
            'parmchk -i mol_charged.mol2 -f mol2 -o mol.frcmod', 'tleap -f leap.in']

    inputs['leap.in'] = '\n'.join(["source leaprc.%s" % ffname,
                                   "%s = loadmol2 mol_charged.mol2" % resname,
                                   "fmod = loadamberparams mol.frcmod",
                                   "check %s" % resname,
                                   "saveoff %s mol.lib" % resname,
                                   "saveamberparm %s mol.prmtop mol.inpcrd" % resname,
                                   "quit\n"])

    def finish_job(j):
        param = ExtraAmberParameters(j.get_output('mol.lib'),
                                     j.get_output('mol.frcmod'),
                                     j)
        return param

    job = pyccc.Job(image=mdt.compute.get_image_path(IMAGE),
                    command=' && '.join(cmds),
                    inputs=inputs,
                    when_finished=finish_job,
                    name="GAFF assignment: %s" % mol.name)

    return mdt.compute.run_job(job, _return_result=True, **kwargs)


