#!/usr/bin/env python
"""
This script drives an NWChem calculation given a generic QM specification
"""
import json
import os
import pint

ureg = pint.UnitRegistry()
ureg.define('bohr = a0 = hbar/(m_e * c * fine_structure_constant')


def convert(q, units):
    """
    Convert a javascript quantity description into a floating point number in the desired units

    Args:
        q (dict): Quantity to convert (of form ``{value: <float>, units: <str>}`` )
        units (str): name of the units

    Returns:
        float: value of the quantity in the passed unit system

    Raises:
        pint.DimensionalityError: If the units are incompatible with the desired quantity

    Examples:
        >>> q = {'value':1.0, 'units':'nm'}
        >>> convert(q, 'angstrom')
        10.0
    """
    quantity = q['value']*ureg(q['units'])
    return quantity.m_as(units)


class QMRunner(object):

    def __init__(self, parameters):
        self.parameters = parameters

    def run_calculation(self):
        """ Drive the calculation, based on passed parameters

        Args:
            parameters (dict): dictionary describing this run (see
               https://github.com/Autodesk/molecular-design-toolkit/wiki/Generic-parameter-names )
        """
        cmd = 'nwchem nw.in > nw.out'
        if 'num_processors' in self.parameters:
            cmd += 'mpirun -n %d' % self.parameters['num_processors']

        os.system(cmd)

    def write_inputs(self):
        """ Write input files using passed parameters

        Args:
            parameters (dict): dictionary describing this run (see
               https://github.com/Autodesk/molecular-design-toolkit/wiki/Generic-parameter-names )
        """
        os.mkdir('./perm')
        inputs = self._make_input_files()
        for filename, contents in inputs.iteritems():
            with open(filename, 'w') as inputfile:
                inputfile.write(contents)


    ##### helper routines below ######
    def _make_input_files(self):
        nwin = [self._header(),
                self._geom_block(),
                self._basisblock(),
                self._chargeblock(),
                self._scfblock(),
                self._otherblocks(),
                self._taskblock()]

        return {'nw.in': '\n'.join(nwin)}

    def _header(self):
        return '\nstart mol\n\npermanent_dir ./perm\nprint medium\n'

    def _geom_block(self):
        if not os.path.isfile('input.xyz'):
            raise IOError('Expecting input coordinates at input.xyz')
        lines = ['geometry units angstrom noautoz noautosym nocenter',
                 '  load format xyz input.xyz',
                 'end']
        return '\n'.join(lines)

    def _basisblock(self):
        return 'basis\n  * library %s\nend' % self.parameters['basis']  # TODO: translate names

    def _scfblock(self):
        lines = [self.SCFTYPES[self.parameters['theory']],
                 self._multiplicityline(),
                 self._theorylines(),
                 'end'
                 ]
        return '\n'.join(lines)

    def _otherblocks(self):
        lines = []
        if self.parameters['runType'] == 'minimization':
            lines.append('driver\n xyz opt\n print high\n')
            if 'minimization_steps' in self.parameters:
                lines.append('maxiter %d' % self.parameters['minimization_steps'])
            lines.append('end\n')

        if 'esp' in self.parameters['properties']:
            lines.append(self._espblock())

        return '\n'.join(lines)

    SCFTYPES = {'rhf': 'scf',
                'uhf': 'scf',
                'rks': 'dft',
                'uks': 'dft',
                'mp2': 'scf'}  # TODO: Not sure if true!!! How do we use DFT reference?

    SCFNAMES = {'rhf': 'rhf',
                'uhf': 'uhf',
                'rks': 'rhf',
                'uks': 'uhf',
                'mp2': 'rhf'}

    TASKNAMES = {'rhf': 'scf',
                 'uhf': 'scf',
                 'rks': 'dft',
                 'uks': 'dft',
                 'mp2': 'mp2'}

    def _scfname(self):
        return self.SCFNAMES[self.parameters['theory']]

    def _taskblock(self):
        if self.parameters['runType'] == 'minimization':
            tasktype = 'optimize'
        elif 'forces' in self.parameters['properties']:
            tasktype = 'gradient'
        else:
            tasktype = 'energy'

        return 'task %s %s' % (self.TASKNAMES[self.parameters['theory']], tasktype)

    STATENAMES = {1: "SINGLET",
                  2: "DOUBLET",
                  3: "TRIPLET",
                  4: "QUARTET",
                  5: "QUINTET",
                  6: "SEXTET",
                  7: "SEPTET",
                  8: "OCTET"}

    def _multiplicityline(self):
        if self.parameters['theory'] in ('rks', 'uks'):
            return 'mult %s' % self.parameters.get('multiplicity', 1)
        else:
            return self.STATENAMES[self.parameters.get('multiplicity', 1)]

    def _constraintblock(self):
        """
        Constraints / restraints are specified in JSON objects at calc.constraints and calc.restraints.

        "Constraints" describe a specific degree of freedom and its value. This value is expected to be
        rigorously conserved by any dynamics or optimization algorithms.

        A "restraint", in constrast, is a harmonic restoring force applied to a degree of freedom.
        Restraint forces should be added to the nuclear hamiltonian. A restraint's "value" is the
        spring's equilibrium position.

        The list at calc['constraints'] has the form:
        [{type: <str>,
          restraint: <bool>,
          value: {value: <float>, units: <str>}
          atomIdx0: <int>...} ...]

        For example, fixed atom constraints have the form:
        {type: 'atomPosition',
         restraint: False,
         atomIdx: <int>}

        Restraints are described similarly, but include a spring constant (of the appropriate units):
        {type: 'bond',
         restraint: True,
         atomIdx1: <int>,
         atomIdx2: <int>,
         springConstant: {value: <float>, units: <str>},
         value: {value: <float>, units: <str>}}
        """
        clist = self.parameters.get('constraints', None)
        if clist is None: return

        lines = ['constraints']

        for constraint in clist:
            if constraint['type'] == 'position':
                lines.append('  fix atom %d' % constraint['atomIdx'])
            elif constraint['type'] == 'distance' and constraint['restraint']:
                k = convert(constraint['springConstant'], 'hartree/(a0*a0)')
                d0 = convert(constraint('value', 'a0'))
                lines.append('  spring bond  %d %d %20.10f %20.10f' %
                             (constraint['atomIdx1']+1, constraint['atomIdx2']+1, k, d0))
            else:
                raise NotImplementedError('Constraint type %s (as restraint: %b) not implemented' %
                                          (constraint['type'], constraint['restraint']))

        lines.append('end')

        return '\n'.join(lines)

    def _espblock(self):
        return "esp\nrecalculate\nend"

    def _chargeblock(self):
        return '\ncharge %s\n' % self.parameters.get('charge', 0)

    def _isdft(self):
        if self.parameters['theory'] in ('rks', 'uks'):
            return True
        elif self.parameters.get('reference', None) in ('rks', 'uks'):
            return True
        else:
            return False

    def _theorylines(self):
        lines = ['direct']
        if self._isdft():
            lines.append('XC %s' % self.parameters['functional'])
        return '\n'.join(lines)


def main(cls):
    with open('params.json', 'r') as pjson:
        parameters = json.load(pjson)
    runner = cls(parameters)
    runner.write_inputs()
    runner.run_calculation()


if __name__ == '__main__':
    main(QMRunner)
