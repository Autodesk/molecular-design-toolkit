#!/usr/bin/env python
"""
This script drives an NWChem calculation given a generic QM specification
"""
import json
import os
import pint

from run import *

CRDPARAMS = "nwchem.crdparams"


class QMMMRunner(QMRunner):
    def _geom_block(self):
        if not os.path.isfile(CRDPARAMS):
            raise IOError('Required input file %s not found' % CRDPARAMS)
        return "mm\ncrdparms load %s\nend" % CRDPARAMS

    def _taskblock(self):
        if self.parameters['runType'] != 'minimization':
            raise ValueError('QMMM only currently supports optimizations')

        if self.parameters['runType'] == 'minimization':
            tasktype = 'optimize'
        elif 'forces' in self.parameters['properties']:
            tasktype = 'gradient'
        else:
            tasktype = 'energy'

        return 'task mm %s %s ignore' % (self._taskname(), tasktype)

if __name__ == '__main__':
    main(QMMMRunner)
