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
