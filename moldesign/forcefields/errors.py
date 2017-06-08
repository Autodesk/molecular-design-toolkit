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

"""
Class that define displays for common errors in assigning a forcefield
"""
import cgi
from .. import utils
from .. import units as u


def show_parameterization_results(errormessages, molin, molout=None):
    print('Forcefield assignment: %s' % ('Success' if molout is not None else 'Failure'))
    for err in errormessages:
        print(utils.html_to_text(err.desc))


class ForcefieldAssignmentError(Exception):
    def __init__(self, messages, mol, molout=None):
        self.errors = messages

    def draw(self):
        # TODO: get this working again
        raise NotImplementedError()



class ForceFieldMessage(object):
    def __init__(self, message):
        self.message = message
        self.desc = message
        self.short = message

    def show(self, viewer):
        viewer.ribbon(opacity=0.7)

    def unshow(self, viewer):
        pass


class StructureOk(ForceFieldMessage):
    """
    A blank message if no other warnings are generated.
    """
    MSG = 'No errors or warnings'

    def __init__(self):
        super().__init__(self.MSG)


class MissingTerms(ForceFieldMessage):
    def __init__(self, message):
        self.message = message
        self.desc = self.message
        self.short = self.message


class UnknownAtom(ForceFieldMessage):
    def __init__(self, message, residue, atom):
        self.message = message
        self.residue = residue
        self.atom = atom
        self.desc = ('ERROR: Atom name "%s" was not found in the "%s" template<br>%s' % (
            self.atom.name, self.residue.resname, self.atom)
                     + '<p>TLeap message:<i>%s</i>' % self.message)
        self.short = 'ERR: %s: unknown atom name "%s" for residue %s' % (self.atom,
                                                                         self.atom.name,
                                                                         self.atom.residue.resname)

    def show(self, viewer):
        viewer.licorice(atoms=self.residue.atoms)
        viewer.vdw(atoms=[self.atom])

    def unshow(self, viewer):
        viewer.ribbon(opacity=0.7)


class MissingAtom(ForceFieldMessage):
    def __init__(self, message, residue, atom):
        self.message = message
        self.residue = residue
        self.atom = atom
        self.desc = 'INFO: Atom %s in %s (chain %s) was added to the system using the "%s" template' % (
            self.atom.name, self.residue.name, self.residue.chain.name, self.residue.resname) + \
                    '<p>TLeap message:<i>%s</i>' % self.message
        self.short = 'INFO: Missing heavy atom %s (index %d)' % (self.atom, self.atom.index)

    def show(self, viewer):
        viewer.licorice(atoms=self.residue.atoms, render=False)
        viewer.vdw(atoms=[self.atom])

    def unshow(self, viewer):
        viewer.ribbon(opacity=0.7)


class UnknownResidue(ForceFieldMessage):
    def __init__(self, message, residue):
        self.message = message
        self.residue = residue
        self._label = None
        self.desc = ('ERROR: Residue type "%s" was not found in residue templates<br>%s' % (self.residue.resname, self.residue)
            + '<p>TLeap message:<i>%s</i>' % self.message)
        self.short = 'ERR: %s: unknown res type "%s"' % (self.residue, self.residue.resname)

    def show(self, viewer):
        viewer.licorice(opacity=1.0, atoms=self.residue.atoms, render=False)
        self._label = viewer.draw_label(position=self.residue.com, text=self.residue.name)

    def unshow(self, viewer):
        viewer.ribbon(opacity=0.7, render=False)
        if self._label: viewer.remove(self._label)
        self._label = None


class UnusualBond(ForceFieldMessage):
    def __init__(self, message, atoms, residues):
        self.message = message
        self.atoms = atoms
        self.residues = residues
        self._has_atoms = (self.atoms[0] is not None) and (self.atoms[1] is not None)
        self._shape = None
        if self._has_atoms:
            self.desc = 'WARNING: Unusual distance between {a1} and {a2}: {d:.3f}'.format(
                d=self.atoms[0].distance(self.atoms[1]), a1=self.atoms[0], a2=self.atoms[1]) \
                        + '<p>TLeap message:<br><i>%s</i>' % cgi.escape(self.message)

            self.short = 'WARN: Unusual dist: {} - {} = ({:.1f})'.format(self.atoms[0],
                                                                   self.atoms[1],
                                                                   self.atoms[0].distance(self.atoms[1]))
        else:
            self.short = 'WARN: Unusual bond - atoms not shown'
            self.desc = 'TLeap message:<br><i>%s</i>' % cgi.escape(self.message)

    def show(self, viewer):
        if self._has_atoms:
            res_opacity = 0.7
        else:
            res_opacity = 0.9

        viewer.licorice(opacity=res_opacity, atoms=self.residues[0], render=False)
        viewer.licorice(opacity=res_opacity, atoms=self.residues[1], render=False)
        if self._has_atoms:
            self._shape = viewer.draw_cylinder(start=self.atoms[0].position,
                                               end=self.atoms[1].position,
                                               radius=0.1 * u.angstrom,
                                               opacity=1.0,
                                               color='red')

    def unshow(self, viewer):
        viewer.ribbon(opacity=0.7, render=False)
        if self._shape: viewer.remove(self._shape)
        self._shape = None