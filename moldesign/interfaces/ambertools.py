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
import moldesign as mdt
from .. import compute, utils
from ..compute import packages
from .. import units as u
from ..molecules import AtomicProperties


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
        charges = {}

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

        mol.properties[chargename] = AtomicProperties(charges)
        return charges

    job = packages.antechamber.make_job(command=command,
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
            mol = mdt.Molecule(mol)
        mdt.helpers.assign_biopolymer_bonds(mol)

        mol.name = '%s-DNA Helix: %s' % (helix_type.upper(), sequence)
        return mol

    job = packages.nab.make_job(command='nab -o buildbdna build.nab && ./buildbdna',
                                inputs={'build.nab': '\n'.join(infile)},
                                name='NAB_build_dna',
                                when_finished=finish_job)

    return mdt.compute.run_job(job, _return_result=True, **kwargs)

