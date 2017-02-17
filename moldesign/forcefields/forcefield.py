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
import moldesign as mdt
from moldesign import units as u


class ForceField(object):
    def __init__(self, mol, amber_params=None, parmed_obj=None):
        if sum(x is not None
               for x in (amber_params, parmed_obj)) > 1:
            raise ValueError('Multiple forcefield definitions passed to %s object'
                             % self.__class__.__name__)

        self.mol = mol
        self.amber_params = self.parmed_obj = self._source_of_truth = None

        if amber_params is not None:
            self.amber_params = amber_params
            self._source_of_truth = 'amber_params'
        elif parmed_obj is not None:
            self.parmed_obj = parmed_obj
            self._source_of_truth = 'parmed_obj'

    def to_parmed(self):
        """ Convert parameters to a parmed structure

        Returns:
            parmed.Structure: a ParmEd object with this parameter set
        """
        if self._source_of_truth == 'parmed_obj':
            return self.parmed_obj
        elif self._source_of_truth == 'amber_params':
            import parmed, os, tempfile

            prmtoppath = os.path.join(tempfile.mkdtemp(), 'prmtop')
            self.amber_params.prmtop.put(prmtoppath)
            pmd = parmed.load_file(prmtoppath,
                                   xyz=self.mol.positions.value_in(u.angstrom))
            return pmd

class FFTerm(object):
    pass


class FFParameters(object):
    """
    This object contains assigned force field parameters for a specific system
    The base focuses on the AMBER / CHARMM - type force
    """

    # TODO: this needs to describe attenuation for things close together
    # TODO: deal with nonbonded exceptions
    def __init__(self, bonds, angles, dihedrals, partial_charges, lennard_jones):
        self.bonds = bonds
        self.angles = angles
        self.dihedrals = dihedrals

        self.partial_charges = partial_charges  # maps atoms to their partial charges
        self.lennard_jones = lennard_jones  # maps atoms to LJ terms

        # lookups
        self.bond_term = {term.bond: term for term in self.bonds}
        self.angle_term = {tuple(term.atoms): term for term in self.angles}
        for term in self.angles: self.angle_term[tuple(reversed(term.atoms))] = term
        self.dihedral_term = {}
        for term in self.dihedrals:
            self.dihedral_term.setdefault(tuple(term.atoms), []).append(term)
            self.dihedral_term.setdefault(tuple(reversed(term.atoms)), []).append(term)


class ResidueTemplate(object):
    def __init__(self, mol, charges, ffparams=None, minimized_mol=None):
        self.mol = mol
        self.charges = {atom.pdbname:charge for atom,charge in charges.iteritems()}

        if ffparams is None:
            self.ffparams = self.get_ffparams()
        else:
            self.ffparams = ffparams

        self.minimized_mol = minimized_mol

    def get_ffparams(self):
        chargemap = {atom: self.charges[atom.pdbname] for atom in self.mol.atoms}
        return mdt.interfaces.ambertools.get_gaff_parameters(self.mol, chargemap)

ffdefaults = dict(protein='ff14SB',
                  dna='OL15',
                  rna='OL3',
                  carbohydrate='GLYCAM_06j-1',
                  lipid='lipid14',
                  water='tip3p',
                  organic='gaff2')
