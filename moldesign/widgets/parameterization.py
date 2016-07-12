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
""" This module contains methods to display 3D information about the assignment of force field
parameters to a biomolecule.

This helps users to visualize errors that occur when trying to assign a forcefield,
such as unrecognized residues, missing atoms, etc.

NOTE:
    This is currently tied to ambertools and tleap! It will need to be made generic if/when
    another method for assigning forcefields is added.
"""
import cgi
import collections
import re

import ipywidgets as ipy

from moldesign import units as u


class ParameterizationDisplay(ipy.Box):
    ATOMSPEC = re.compile(r'\.R<(\S+) (\d+)>\.A<(\S+) (\d+)>')

    def __init__(self, job, molin, molout=None):
        self.molin = molin
        self.molout = molout
        self.job = job
        self.msg = []
        self.parse_errors(self.job)

        self.status = ipy.HTML('<h4>Forcefield assignment: %s</h4>' %
                               ('Success' if molout else 'FAILED'))

        self.listdesc = ipy.HTML('<b>Errors / warnings:</b>')
        self.errorlist = ipy.Select(options=collections.OrderedDict((e.short, e) for e in self.msg))
        self.errmsg = ipy.HTML('-')

        self.viewerpane = self.molin.draw()
        self.viewer = self.viewerpane.children[0]
        self.viewer.ribbon(opacity=0.7)

        if self.errorlist.value is not None:
            self.switch_display({'old': self.errorlist.value, 'new': self.errorlist.value})
        self.errorlist.observe(self.switch_display, 'value')
        children = (self.status,
                    ipy.HBox([self.viewerpane, ipy.VBox([self.listdesc, self.errorlist])]),
                    self.errmsg)

        super(ParameterizationDisplay, self).__init__(children=children)

    def switch_display(self, d):
        old = d['old']
        old.unshow(self.viewer)
        self.errmsg.value = '-'
        new = d['new']
        new.show(self.viewer)
        self.errmsg.value = new.desc

    def parse_errors(self, job):
        # TODO: special messages for known problems (e.g. histidine)
        unknown_res = set()
        lineiter = iter(job.stdout.split('\n'))
        reslookup = {str(i + self.molin.residues[0].pdbindex): r for i,r in enumerate(self.molin.residues)}

        def _atom_from_re(s):
            resname, residx, atomname, atomidx = s
            r = reslookup[residx]
            a = r[atomname]
            return a

        def unusual_bond(l):
            atomre1, atomre2 = self.ATOMSPEC.findall(l)
            try:
                a1, a2 = _atom_from_re(atomre1), _atom_from_re(atomre2)
            except KeyError:
                a1 = a2 = None
            r1 = reslookup[atomre1[1]]
            r2 = reslookup[atomre2[1]]
            self.msg.append(UnusualBond(l, (a1, a2), (r1, r2)))

        while True:
            try: line = lineiter.next()
            except StopIteration: break

            fields = line.split()
            if fields[0:2] == ['Unknown','residue:']:
                # EX: "Unknown residue: 3TE   number: 499   type: Terminal/beginning"
                res = self.molin.residues[int(fields[4])]
                self.msg.append(UnknownResidue(line,res))
                unknown_res.add(res)

            elif fields[:4] == 'Warning: Close contact of'.split():
                # EX: "Warning: Close contact of 1.028366 angstroms between .R<DC5 1>.A<HO5' 1> and .R<DC5 81>.A<P 9>"
                unusual_bond(line)

            elif fields[:6] == 'WARNING: There is a bond of'.split():
                # Matches two lines, EX:
                # "WARNING: There is a bond of 34.397700 angstroms between:"
                # "-------  .R<DG 92>.A<O3' 33> and .R<DG 93>.A<P 1>"
                nextline = lineiter.next()
                unusual_bond(line + nextline)

            elif fields[:5] == 'Created a new atom named:'.split():
                # EX: "Created a new atom named: P within residue: .R<DC5 81>"
                residue = reslookup[fields[-1][:-1]]
                if residue in unknown_res: continue  # suppress atoms from an unknown res ...
                atom = residue[fields[5]]
                self.msg.append(UnknownAtom(line, residue, atom))


class ForceFieldMessage(object):
    pass


class UnknownAtom(ForceFieldMessage):
    def __init__(self, message, residue, atom):
        self.message = message
        self.residue = residue
        self.atom = atom
        self.desc = 'ERROR: Atom name "%s" was not found in the "%s" template<br>%s' % (
            self.atom.name, self.residue.resname, self.atom) + '<p>TLeap message:<i>%s</i>' % self.message
        self.short = 'ERR: %s: unknown atom name "%s" for residue %s' % (self.atom,
                                                                         self.atom.name,
                                                                         self.atom.residue.resname)

    def show(self, viewer):
        viewer.licorice(atoms=self.residue.atoms, render=False)
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
        else:
            viewer.render()

    def unshow(self, viewer):
        viewer.ribbon(opacity=0.7, render=False)
        if self._shape: viewer.remove(self._shape)
        self._shape = None