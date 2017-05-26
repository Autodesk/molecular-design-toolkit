#!/usr/bin/env python
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
""" Generates a residue database for MDT from the PDB Chemical Component Dictionary

Due to its size, we generate a dbm database that can be accessed without loading the whole thing
into memory.
"""
from __future__ import print_function
from builtins import zip
import collections
import itertools
import sys
import os
sys.path.insert(0, '..')

from moldesign import utils

# to install mmCif: pip install git+git://github.com/glenveegee/PDBeCIF
from mmCif.mmcifIO import CifFileReader

ORDERS = dict(SING=1, DOUB=2, TRIP=3)

DBNAME = 'chemical_components'


def bond_data(data):
    if '_chem_comp_bond' not in data:
        return {}

    bonds = collections.OrderedDict()
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
    return bonds

ATOMFIELDS = ['symbol', 'formal_charge']


def atom_data(data):
    atomdata = data['_chem_comp_atom']

    atoms = collections.OrderedDict()
    if isinstance(atomdata['atom_id'], list):
        atom_ids = atomdata['atom_id']
        atom_charges = atomdata['charge']
        atom_symbols = atomdata['type_symbol']
    else:
        atom_ids = [atomdata['atom_id']]
        atom_charges = [atomdata['charge']]
        atom_symbols = [atomdata['type_symbol']]

    for id, charge, symbol in zip(atom_ids, atom_charges, atom_symbols):
        atoms[id] = [symbol, int(charge)]

    return atoms


def download_ccd():
    os.system('wget -N ftp://ftp.wwpdb.org/pub/pdb/data/monomers/aa-variants-v1.cif.gz')
    os.system('wget -N ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif.gz')
    os.system('gunzip -v -f -k components.cif.gz aa-variants-v1.cif.gz')

RESFIELDS = ['name', 'type', 'atoms', 'bonds']


def read_ccd():
    print('Reading components.cif ...')
    sys.stdout.flush()
    components = CifFileReader()
    all_residues = components.read('components.cif')
    print('Reading aa-variants-v1.cif ...')
    amino_variants = components.read('aa-variants-v1.cif')
    return itertools.chain(iter(all_residues.items()), iter(amino_variants.items()))


def store_residue(resname, data, db):
    try:
        record = [data['_chem_comp']['name'],
                  data['_chem_comp']['type'],
                  atom_data(data),
                  bond_data(data)]
    except (KeyError, ValueError) as exc:
        print('\nFAILED for residue %s: %s'%(resname, exc))
    else:
        db[resname] = record
    print(resname, end=' ')
    sys.stdout.flush()


def main():
    print('This program regenerates the `%s` database using the "Chemical Component ' % DBNAME)
    print('Dictionary" and "Protonation Variants Dictionary" from http://www.wwpdb.org/data/ccd')

    if len(sys.argv) > 1 and sys.argv[1] == '--download':
        print('Downloading newest CCD files ...\n')
        download_ccd()
    else:
        print('Using cached download (run with "--download" to download a new version)')

    residues = read_ccd()

    db = utils.CompressedJsonDbm(DBNAME, 'n')
    db['__FIELDS__'] = {'ATOMFIELDS': ATOMFIELDS,
                        'RESFIELDS': RESFIELDS}

    for resname, data in residues:
        store_residue(resname, data, db)

    print('\nCreated db "%s" (%d records)' % (DBNAME, len(db)))
    db.close()


if __name__ == '__main__':
    main()
