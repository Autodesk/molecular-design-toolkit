from __future__ import print_function, absolute_import, division
from future.builtins import *
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
from past.builtins import basestring
import re
import os
import tempfile

import moldesign as mdt
import pyccc
from .. import compute, utils
from .. import units as u
from .. import forcefields

IMAGE = 'ambertools'


class AmberParameters(object):
    def __getstate__(self):
        state = self.__dict__.copy()
        state['job'] = None
        return state

    def __init__(self, prmtop, inpcrd, job=None):
        self.prmtop = prmtop
        self.inpcrd = inpcrd
        self.job = job

    def to_parmed(self):
        import parmed
        prmtoppath = os.path.join(tempfile.mkdtemp(), 'prmtop')
        self.prmtop.put(prmtoppath)
        pmd = parmed.load_file(prmtoppath)
        return pmd


class ExtraAmberParameters(object):
    def __getstate__(self):
        state = self.__dict__.copy()
        state['job'] = None
        return state

    def __init__(self, lib, frcmod, job=None):
        self.lib = lib
        self.frcmod = frcmod
        self.job = job


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
    command = 'antechamber -fi mol2 -i mol.mol2 -fo mol2 -o out.mol2 -c %s -an n'%ambname
    if charge != 0:
        command += ' -nc %d' % charge.value_in(u.q_e)

    def finish_job(job):
        """Callback to complete the job"""
        lines = iter(job.get_output('out.mol2').read().split('\n'))
        charges = utils.DotDict(type='atomic')

        line = next(lines)
        while line.strip()[:len('@<TRIPOS>ATOM')] != '@<TRIPOS>ATOM':
            line = next(lines)

        line = next(lines)
        while line.strip()[:len('@<TRIPOS>BOND')] != '@<TRIPOS>BOND':
            fields = line.split()
            idx = int(fields[0])-1
            assert mol.atoms[idx].name == fields[1]
            charges[mol.atoms[idx]] = u.q_e*float(fields[-1])
            line = next(lines)

        mol.properties[chargename] = charges
        return charges

    job = pyccc.Job(image=mdt.compute.get_image_path(IMAGE),
                    command=command,
                    name="%s, %s" % (chargename, mol.name),
                    inputs={'mol.mol2': mol.write(format='mol2')},
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
                    image=mdt.compute.get_image_path('nucleic_acid_builder'),
                    inputs={'build.nab': '\n'.join(infile)},
                    name='NAB_build_dna',
                    when_finished=finish_job)

    return mdt.compute.run_job(job, _return_result=True, **kwargs)


@utils.kwargs_from(compute.run_job)
def _run_tleap_assignment(mol, leapcmds, files=None, **kwargs):
    """
    Drives tleap to create a prmtop and inpcrd file. Specifically uses the AmberTools 16
    tleap distribution.

    Defaults are as recommended in the ambertools manual.

    Args:
        mol (moldesign.Molecule): Molecule to set up
        leapcmds (List[str]): list of the commands to load the forcefields
        files (List[pyccc.FileReference]): (optional) list of additional files
           to send
        **kwargs: keyword arguments to :meth:`compute.run_job`

    References:
        Ambertools Manual, http://ambermd.org/doc12/Amber16.pdf. See page 33 for forcefield
        recommendations.
    """
    leapstr = leapcmds[:]
    inputs = {}
    if files is not None:
        inputs.update(files)
    inputs['input.pdb'] = mol.write(format='pdb')

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


def _prep_for_tleap(mol):
    """ Returns a modified *copy* that's been modified for input to tleap

    Makes the following modifications:
       1. Reassigns all residue IDs
       2. Assigns tleap-appropriate cysteine resnames
    """
    change = False
    clean = mdt.Molecule(mol.atoms)
    for residue in clean.residues:
        residue.pdbindex = residue.index+1

        if residue.resname == 'CYS':  # deal with cysteine states
            if 'SG' not in residue.atoms or 'HG' in residue.atoms:
                continue  # sulfur's missing, we'll let tleap create it
            else:
                sulfur = residue.atoms['SG']

            if sulfur.formal_charge == -1*u.q_e:
                residue.resname = 'CYM'
                change = True
                continue

            # check for a reasonable hybridization state
            if sulfur.formal_charge != 0 or sulfur.num_bonds not in (1, 2):
                raise ValueError("Unknown sulfur hybridization state for %s"
                                 % sulfur)

            # check for a disulfide bond
            for otheratom in sulfur.bonded_atoms:
                if otheratom.residue is not residue:
                    if otheratom.name != 'SG' or otheratom.residue.resname not in ('CYS', 'CYX'):
                        raise ValueError('Unknown bond from cysteine sulfur (%s)' % sulfur)

                    # if we're here, this is a cystine with a disulfide bond
                    print('INFO: disulfide bond detected. Renaming %s from CYS to CYX' % residue)
                    sulfur.residue.resname = 'CYX'

            clean._rebuild_topology()

    return clean


@utils.kwargs_from(mdt.compute.run_job)
def create_ff_parameters(mol, charges='esp', baseff='gaff2', **kwargs):
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
        baseff (str): Name of the gaff-like forcefield file (default: gaff2)

    Returns:
        TLeapForcefield: Forcefield object for this residue
    """
    # Check that there's only 1 residue, give it a name
    assert mol.num_residues == 1
    if mol.residues[0].resname is None:
        mol.residues[0].resname = 'UNL'
        print('Assigned residue name "UNL" to %s' % mol)
    resname = mol.residues[0].resname

    # check that atoms have unique names
    if len(set(atom.name for atom in mol.atoms)) != mol.num_atoms:
        raise ValueError('This molecule does not have uniquely named atoms, cannot assign FF')

    if charges == 'am1-bcc' and 'am1-bcc' not in mol.properties:
        calc_am1_bcc_charges(mol)
    elif charges == 'gasteiger' and 'gasteiger' not in mol.properties:
        calc_gasteiger_charges(mol)
    elif charges == 'esp' and 'esp' not in mol.properties:
        # TODO: use NWChem ESP to calculate
        raise NotImplementedError()

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

    cmds = ['antechamber -i mol.mol2 -fi mol2 -o mol_charged.mol2 '
                   ' -fo mol2 -c rc -cf mol.charges -rn %s' % resname,
            'parmchk -i mol_charged.mol2 -f mol2 -o mol.frcmod',
            'tleap -f leap.in',
            'sed -e "s/tempresname/%s/g" mol_rename.lib > mol.lib' % resname]

    base_forcefield = forcefields.TLeapLib(baseff)
    inputs['leap.in'] = '\n'.join(["source leaprc.%s"%baseff,
                                   "tempresname = loadmol2 mol_charged.mol2",
                                   "fmod = loadamberparams mol.frcmod",
                                   "check tempresname",
                                   "saveoff tempresname mol_rename.lib",
                                   "saveamberparm tempresname mol.prmtop mol.inpcrd",
                                   "quit\n"])

    def finish_job(j):
        leapcmds = ['source leaprc.gaff2']
        files = {}
        for fname, f in j.glob_output("*.lib").items():
            leapcmds.append('loadoff %s' % fname)
            files[fname] = f
        for fname, f in j.glob_output("*.frcmod").items():
            leapcmds.append('loadAmberParams %s' % fname)
            files[fname] = f

        param = forcefields.TLeapForcefield(leapcmds, files)
        param.add_ff(base_forcefield)
        param.assign(mol)
        return param

    job = pyccc.Job(image=mdt.compute.get_image_path(IMAGE),
                    command=' && '.join(cmds),
                    inputs=inputs,
                    when_finished=finish_job,
                    name="GAFF assignment: %s" % mol.name)

    return mdt.compute.run_job(job, _return_result=True, **kwargs)


ATOMSPEC = re.compile(r'\.R<(\S+) ([\-0-9]+)>\.A<(\S+) ([\-0-9]+)>')


def _parse_tleap_errors(job, molin):
    # TODO: special messages for known problems (e.g. histidine)
    msg = []
    unknown_res = set()  # so we can print only one error per unkonwn residue
    lineiter = iter(job.stdout.split('\n'))
    offset = utils.if_not_none(molin.residues[0].pdbindex, 1)
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
        return forcefields.errors.UnusualBond(l, (a1, a2), (r1, r2))

    def _parse_tleap_logline(line):
        fields = line.split()
        if fields[0:2] == ['Unknown', 'residue:']:
            # EX: "Unknown residue: 3TE   number: 499   type: Terminal/beginning"
            res = molin.residues[int(fields[4])]
            unknown_res.add(res)
            return forcefields.errors.UnknownResidue(line, res)

        elif fields[:4] == 'Warning: Close contact of'.split():
            # EX: "Warning: Close contact of 1.028366 angstroms between .R<DC5 1>.A<HO5' 1> and .R<DC5 81>.A<P 9>"
            return unusual_bond(line)

        elif fields[:6] == 'WARNING: There is a bond of'.split():
            # Matches two lines, EX:
            # "WARNING: There is a bond of 34.397700 angstroms between:"
            # "-------  .R<DG 92>.A<O3' 33> and .R<DG 93>.A<P 1>"
            nextline = next(lineiter)
            return unusual_bond(line+nextline)

        elif fields[:5] == 'Created a new atom named:'.split():
            # EX: "Created a new atom named: P within residue: .R<DC5 81>"
            residue = reslookup[fields[-1][:-1]]
            if residue in unknown_res:
                return None  # suppress atoms from an unknown res ...
            atom = residue[fields[5]]
            return forcefields.errors.UnknownAtom(line, residue, atom)

        elif (fields[:5] == '** No torsion terms for'.split() or
                      fields[:5] == 'Could not find angle parameter:'.split()):
            # EX: " ** No torsion terms for  ca-ce-c3-hc"
            return forcefields.errors.MissingTerms(line.strip())

        else:  # ignore this line
            return None

    while True:
        try:
            line = next(lineiter)
        except StopIteration:
            break

        try:
            errmsg = _parse_tleap_logline(line)
        except (KeyError, ValueError):
            print("WARNING: failed to process TLeap message '%s'" % line)
            msg.append(forcefields.errors.ForceFieldMessage(line))

        else:
            if errmsg is not None:
                msg.append(errmsg)

    return msg


