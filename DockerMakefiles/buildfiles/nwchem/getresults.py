#!/usr/bin/env python
"""
This script extracts data from the NWChem DB and outputs it in a standardized
JSON format. It uses Marat Valiev's simple API from https://github.com/nwchem-python/nwapi

Note:
    This DOES NOT handle general cases! We make some very implicit assumptions about the nature
    of the calculation; the parsed file is likely to be inaccurate if it's not the result of
    an MDT-driven calculation
"""
import json
import os
import nwchem

FORMAT = {'format': 'MolecularDesignInterchange',
          'version': '0.8alpha'}


##### Software metadata ####

def get_package_name(): return 'NWChem6.6'

def get_package_version():
    filefields = os.environ['NWCHEMFILE'].split()

    while filefields[-1] in 'gz zip tar bz2 xz'.split():
        filefields.pop()

    return '.'.join(filefields)

def get_package_citation():
    return 'doi:10.1016/j.cpc.2010.04.018'

def get_docker_image(): # TODO: do this properly
    return 'docker.io/autodesk/moldesign:nwchem-[MDTVERSION]'

def prep():
    global _PROPERTYGROUP

    nwchem.rtdb_open('perm/mol.db', 'old')

    _PROPERTYGROUP = nwchem.rtdb_get('task:theory')


##### Units #####

def get_position_units():
    return 'a0'

def get_dipole_units():
    return 'q_e*a0'

def get_mass_units():
    return 'amu'

def get_energy_units():
    return 'hartree'

def get_force_units():
    return 'hartree/a0'


##### Calculation description #####

# TODO: semantics for reference wfns / post-HF methods. How do we denote CIS/UHF?

def get_theory():
    theory = nwchem.rtdb_get('task:theory').lower()
    if theory == 'dft':
        scf = nwchem.rtdb_get('dft:scftype').lower()
    else:
        scf = nwchem.rtdb_get('scf:scftype').lower()

    if theory == 'scf':
        return scf
    elif theory == 'dft':
        if scf == 'uhf':
            return 'uks'
        elif scf =='rhf':
            return 'rks'

    # if here, it wasn't handled
    raise NotImplementedError('%s %s' % (theory, scf))


def get_scftype():
    if nwchem.rtdb_get('task:theory').lower() == 'dft':
        return nwchem.rtdb_get('dft:scftype')
    else:
        return nwchem.rtdb_get('scf:scftype')

def get_basis():
    return nwchem.rtdb_get('basis:ao basis:star bas type')


##### Molecular information #####

def get_symmetry():
    return nwchem.rtdb_get('geometry:geometry:group name')

def get_spin_multiplicity():
    return nwchem.rtdb_get('dft:mult')

def get_total_charge():
    return int(nwchem.rtdb_get('charge'))

def get_positions():
    return _reshape_atom_vector(nwchem.rtdb_get('geometry:geometry:coords'))

def get_atomic_numbers():
    return map(int, nwchem.rtdb_get('geometry:geometry:charges'))

def get_atom_masses():
    return nwchem.rtdb_get('geometry:geometry:masses')

def get_atom_names():
    return nwchem.rtdb_get('geometry:geometry:tags')

def get_functional():
    return nwchem.rtdb_get('dft:xc_spec')


##### Caculated quantities #####

# TODO: need to parse ouput file for electronic energy, nuclear repulsion, 1-electron energy, etc.
# TODO: orbitals, overlaps, fock matrix, ... (.movecs, .aoints files)

def get_converged():
    return 1 == nwchem.rtdb_get('%s:converged' % _PROPERTYGROUP)

def get_dipole():
    return nwchem.rtdb_get('%s:dipole' % _PROPERTYGROUP)

def get_forces():
    try:
        forces = nwchem.rtdb_get('%s:gradient' % _PROPERTYGROUP)
        return _reshape_atom_vector([-f for f in forces])
    except:
        return None

def get_potential_energy():
    return nwchem.rtdb_get('%s:energy' % _PROPERTYGROUP)


# geometry:driverinitial:coords - these are the initial coordinates for optimization, presumably


################################# helper functions below ###########################
def _reshape_atom_vector(v):
    assert len(v) % 3 == 0
    return [v[i:i+3] for i in xrange(0,len(v), 3)]


def _get_method_description():
    result = {
        'package': {'name': get_package_name(),
                    'version': get_package_version(),
                    'citation': get_package_citation()},
        'basis': get_basis(),
        'scf': get_scftype(),
        'theory': get_theory()}

    if result['theory'] == 'dft':
        result['functional'] = get_functional()

    return result

def _get_calculated_properties():
    props = {'method': _get_method_description()}

    _insert_if_present(props, 'dipole', get_dipole(), get_dipole_units())
    _insert_if_present(props, 'forces', get_forces(), get_force_units())
    _insert_if_present(props, 'potential_energy', get_potential_energy(), get_energy_units())

    return props


def _insert_if_present(d, key, val, unitval=None):
    if val is not None:
        if unitval is None:
            d[key] = val
        else:
            d[key] = {'value': val, 'units': unitval}


def _main():
    prep()

    state = {'positions': {'value': get_positions(),
                           'units': get_position_units()},
             'calculated': _get_calculated_properties(),

             'charge': get_total_charge(),
             'spin_multiplicity': get_spin_multiplicity(),
             'symmetry': get_symmetry()
             }

    topology = {'atomArray':
                    {'names': get_atom_names(),
                     'masses': {'value': get_atom_masses(),
                                'units': get_mass_units()},
                     'atomicNumbers': get_atomic_numbers()}
                }

    result = {'format': FORMAT,
              'topology': topology,
              'states': [state]}

    with open('results.json', 'w') as outfile:
        json.dump(result, outfile)


if __name__ == '__main__':
    _main()
