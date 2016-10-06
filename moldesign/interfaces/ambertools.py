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
import re

import traitlets

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
    print('DeprecationWarning: build_bdna is deprecated. '
          "Use `build_dna_helix(sequence, helix_type='b')` instead")
    return build_dna_helix(sequence, helix_type='b', **kwargs)


@utils.kwargs_from(mdt.compute.run_job)
def build_dna_helix(sequence, helix_type='B', **kwargs):
    """ Uses Ambertools' Nucleic Acid Builder to build a 3D DNA double-helix.

    Args:
        sequence (str): DNA sequence for one of the strands (a complementary sequence will
            automatically be created)
        helix_type (str): Type of helix - 'A'=Arnott A-DNA
                                          'B'=B-DNA (from standard templates and helical params),
                                          'LB'=Langridge B-DNA,
                                          'AB'=Arnott B-DNA,
                                          'SB'=Sasisekharan left-handed B-DNA
        **kwargs: arguments for :meth:`compute.run_job`

    All helix types except 'B' are taken from fiber diffraction data (see the refernce for details)

    Returns:
        moldesign.Molecule: B-DNA double helix

    References:
        See NAB / AmberTools documentation: http://ambermd.org/doc12/Amber16.pdf, pg 771-2
    """
    infile = ['molecule m;']
    if helix_type.lower() == 'b':
        infile.append('m = bdna( "%s" );' % sequence.lower())
    else:
        infile.append('m = fd_helix( "%sdna",  "%s", "dna" );'
                      % (helix_type.lower(), sequence.lower()))
    infile.append('putpdb( "helix.pdb", m, "-wwpdb");\n')

    def finish_job(job):
        mol = mdt.fileio.read_pdb(job.get_output('helix.pdb').open(), assign_ccd_bonds=False)
        if mol.num_chains == 1:
            assert mol.num_residues % 2 == 0
            oldchain = mol.chains[0]
            oldchain.name = oldchain.pdbindex = oldchain.pdbname = 'A'
            newchain = mdt.Chain('B')
            for residue in mol.residues[mol.num_residues//2:]:
                residue.chain = newchain
                for atom in residue:
                    atom.chain = newchain
            mol = mdt.Molecule(mol)
        mdt.helpers.assign_biopolymer_bonds(mol)

        mol.name = '%s-DNA Helix: %s' % (helix_type.upper(), sequence)
        return mol

    job = pyccc.Job(command='nab -o buildbdna build.nab && ./buildbdna',
                    image=mdt.compute.get_image_path(IMAGE),
                    inputs={'build.nab': '\n'.join(infile)},
                    name='NAB_build_dna',
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

    errors = _parse_tleap_errors(job, mol)

    try:
        report = ParameterizationDisplay(errors, mol, molout=newmol)
        uibase.display_log(report, title='ERRORS/WARNINGS', show=True)
    except traitlets.TraitError:
        print 'Forcefield assignment: %s' % ('Success' if newmol is not None else 'Failure')
        for err in errors:
            print err.desc

    if newmol is not None:
        return newmol
    else:
        raise ParameterizationError('TLeap failed to assign force field parameters for %s' % mol, job)


@utils.kwargs_from(mdt.compute.run_job)
def parameterize(mol, charges='esp', ffname='gaff2', **kwargs):
    """Parameterize ``mol``, typically using GAFF parameters.

    This will both assign a forcefield to the molecule (at ``mol.ff``) and produce the parameters
    so that they can be used in other systems (e.g., so that this molecule can be simulated
    embedded in a larger protein)

    Note:
        'am1-bcc' and 'gasteiger' partial charges will be automatically computed if necessary.
        Other charge types must be precomputed.

    Args:
        mol (moldesign.Molecule):
        charges (str or dict): what partial charges to use? Can be a dict (``{atom:charge}``) OR
            a string, in which case charges will be read from
           ``mol.properties.[charges name]``; typical values will be 'esp', 'mulliken',
           'am1-bcc', etc. Use 'zero' to set all charges to 0 (for QM/MM and testing)
        ffname (str): Name of the gaff-like forcefield file (default: gaff2)

    Returns:
        ExtraAmberParameters: Parameters for the molecule; this object can be used to create
            forcefield parameters for other systems that contain this molecule
    """
    assert mol.num_residues == 1
    if mol.residues[0].resname is None:
        mol.residues[0].resname = 'UNL'
        print 'Assigned residue name "UNL" to %s' % mol
    resname = mol.residues[0].resname

    if charges == 'am1-bcc' and 'am1-bcc' not in mol.properties:
        calc_am1_bcc_charges(mol)
    elif charges == 'gasteiger' and 'gasteiger' not in mol.properties:
        calc_gasteiger_charges(mol)

    if charges == 'zero':
        charge_array = [0.0 for atom in mol.atoms]
    elif isinstance(charges, basestring):
        charge_array = u.array([mol.properties[charges][atom] for atom in mol.atoms])
        if not charge_array.dimensionless:  # implicitly convert floats to fundamental charge units
            charge_array = charge_array.to(u.q_e).magnitude
    else:
        charge_array = [charges[atom] for atom in mol.atoms]

    inputs = {'mol.mol2': mol.write(format='mol2'),
              'mol.charges': '\n'.join(map(str, charge_array))}

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
        tempmol = mdt.assign_forcefield(mol, parameters=param)
        mol.ff = tempmol.ff
        return param

    job = pyccc.Job(image=mdt.compute.get_image_path(IMAGE),
                    command=' && '.join(cmds),
                    inputs=inputs,
                    when_finished=finish_job,
                    name="GAFF assignment: %s" % mol.name)

    return mdt.compute.run_job(job, _return_result=True, **kwargs)


ATOMSPEC = re.compile(r'\.R<(\S+) (\d+)>\.A<(\S+) (\d+)>')


def _parse_tleap_errors(job, molin):
    from moldesign.widgets.parameterization import UnusualBond, UnknownAtom, UnknownResidue

    # TODO: special messages for known problems (e.g. histidine)
    msg = []
    unknown_res = set()  # so we can print only one error per unkonwn residue
    lineiter = iter(job.stdout.split('\n'))
    offset = utils.if_not_none(molin.residues[0].pdbindex, 0)
    reslookup = {str(i+offset): r for i,r in enumerate(molin.residues)}

    def _atom_from_re(s):
        resname, residx, atomname, atomidx = s
        r = reslookup[residx]
        a = r[atomname]
        return a

    def unusual_bond(l):
        atomre1, atomre2 = ATOMSPEC.findall(l)
        try:
            a1, a2 = _atom_from_re(atomre1), _atom_from_re(atomre2)
        except KeyError:
            a1 = a2 = None
        r1 = reslookup[atomre1[1]]
        r2 = reslookup[atomre2[1]]
        return UnusualBond(l, (a1, a2), (r1, r2))

    while True:
        try:
            line = lineiter.next()
        except StopIteration:
            break

        fields = line.split()
        if fields[0:2] == ['Unknown','residue:']:
            # EX: "Unknown residue: 3TE   number: 499   type: Terminal/beginning"
            res = molin.residues[int(fields[4])]
            unknown_res.add(res)
            msg.append(UnknownResidue(line,res))

        elif fields[:4] == 'Warning: Close contact of'.split():
            # EX: "Warning: Close contact of 1.028366 angstroms between .R<DC5 1>.A<HO5' 1> and .R<DC5 81>.A<P 9>"
            msg.append(unusual_bond(line))

        elif fields[:6] == 'WARNING: There is a bond of'.split():
            # Matches two lines, EX:
            # "WARNING: There is a bond of 34.397700 angstroms between:"
            # "-------  .R<DG 92>.A<O3' 33> and .R<DG 93>.A<P 1>"
            nextline = lineiter.next()
            msg.append(unusual_bond(line + nextline))

        elif fields[:5] == 'Created a new atom named:'.split():
            # EX: "Created a new atom named: P within residue: .R<DC5 81>"
            residue = reslookup[fields[-1][:-1]]
            if residue in unknown_res: continue  # suppress atoms from an unknown res ...
            atom = residue[fields[5]]
            msg.append(UnknownAtom(line, residue, atom))


    return msg