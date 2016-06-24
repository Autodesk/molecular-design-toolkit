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
#!/usr/bin/env python
""" Script to generate connection tables from the official PDB definitions

Due to its size, we generate a dbm database that can be accessed without loading the whole thing
into memory.
"""
import itertools
import sys
sys.path.insert(0, '..')
import utils

# to install mmCif: pip install git+git://github.com/glenveegee/PDBeCIF
from mmCif.mmcifIO import CifFileReader

ORDERS = dict(SING=1, DOUB=2, TRIP=3)


if __name__ == '__main__':
    print 'This program regenerates the `residue_bonds` database using the "Chemical Component '
    print 'Dictionary" and "Protonation Variants Dictionary" from http://www.wwpdb.org/data/ccd'
    print 'These files are not included in the repository and must be downloaded manually.\n'
    print 'Reading components.cif ...'
    sys.stdout.flush()

    components = CifFileReader()
    all_residues = components.read('components.cif')

    print 'Reading aa-variants-v1.cif ...'
    amino_variants = components.read('aa-variants-v1.cif')

    residues = []

    db = utils.CompressedJsonDbm('residue_bonds', 'n')

    for resname, data in itertools.chain(all_residues.iteritems(), amino_variants.iteritems()):
        if '_chem_comp_bond' not in data:
            print '\nskipped %s (no bonds)' % resname
            continue

        bonds = {}
        bond_data = data['_chem_comp_bond']
        if type(bond_data['atom_id_1']) == str:  # there's just one bond
            a1 = bond_data['atom_id_1']
            a2 = bond_data['atom_id_2']
            order = ORDERS[bond_data['value_order']]
            bonds[a1] = {a2: order}
            bonds[a2] = {a1: order}
        else:
            for a1, a2, order in zip(bond_data['atom_id_1'],
                                     bond_data['atom_id_2'],
                                     bond_data['value_order']):
                bonds.setdefault(a1, {})[a2] = ORDERS[order]
                bonds.setdefault(a2, {})[a1] = ORDERS[order]

        db[resname] = bonds

        print resname,
        sys.stdout.flush()

    print 'Created db (%d records)' % len(db)
    db.close()







