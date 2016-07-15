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

try:
    import pdbfixer as pf
except ImportError:
    force_remote = True
    pf = None
else:
    force_remote = False

from moldesign import units as u

from .  import openmm as opm


def add_hydrogen(mol, pH=7.4):
    fixer = _get_fixer(mol)
    fixer.addMissingHydrogens(pH)
    return _get_mol(fixer)


def mutate(mol, mutationmap):
    fixer = _get_fixer(mol)
    chain_mutations = {}
    for res in mutationmap:
        chain_mutations.setdefault(res.chain.pdbname, {})[res] = mutationmap[res]

    for chainid, mutations in chain_mutations.iteritems():
        mutstrings = ['%s-%d-%s' % (res.resname, res.pdbindex, newname) for res,newname in mutations.iteritems()]
        print 'Applying mutations to chain %s: %s' % (chainid, ', '.join(mutstrings))
        fixer.applyMutations(mutations, chainid)
    return _get_mol(fixer)


def solvate(mol, padding=10.0*u.angstrom, size=None, cation='Na+', anion='Cl-'):
    fixer = _get_fixer(mol)
    if size is not None: size = size.to_simtk()
    if padding is not None: padding = padding.to_simtk()
    fixer.addSolvent(padding=padding, boxSize=size, positiveIon=cation, negativeIon=anion)
    return _get_mol(fixer)


def _get_fixer(mol):
    mol.write('/tmp/tmp.pdb', format='pdb')
    fixer = pf.PDBFixer('/tmp/tmp.pdb')
    return fixer

def _get_mol(f):
    return opm.topology_to_mol(f.topology, positions=f.positions)

