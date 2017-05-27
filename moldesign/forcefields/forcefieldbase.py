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

import moldesign as mdt
from .errors import ForcefieldAssignmentError, show_parameterization_results
from ..utils import exports


class Forcefield(object):
    """ Abstract base class representing a biomolecular forcefield definition,
    such as amber14sb or charm22.

    These contain both atom type templates AND forcefield parameters.
    """

    def assign(self, mol, display=True):
        """ Assign this forcefield to a molecule.

        This will create a new ForcefieldParams object at ``mol.ff``. The molecule will not be
        otherwise modified

        Args:
            mol (moldesign.Molecule): molecule to assign the forcefield to. It's ``ff``
               attribute will be overwritten with the new parameters
            display (bool): in jupyter, creates a parameterization display if the assignment fails
               (ignored outside of jupyter)

        Raises:
            ForcefieldAssignmentError: if the forcefield cannot be assigned
        """
        raise NotImplementedError()

    def create_prepped_molecule(self, mol, display=True):
        """ Return a modified _copy_ of the passed molecule with assigned forcefield.

        The modified copy may be modified in one or more of the following ways:
          1) Missing atoms added from residue templates
          2) Bond topology changed to fit residue templates
          3) Residue names changed to fit this FF's residue naming conventions

        Args:
            mol (moldesign.Molecule): molecule to create a copy of for forcefield assignment
            display (bool): in jupyter, creates a parameterization error display upon failure
               (ignored outside of jupyter)

        Raises:
            ForcefieldAssignmentError: if the forcefield cannot be assigned

        Returns:
            moldesign.Molecule: A modified copy of the passed molecule
        """
        raise NotImplementedError()

    def check(self, mol, display=True):
        """ Return True if this forcefield can be assigned to the given molecule without errors,
        otherwise return the exception object (without raising it)

        Args:
            mol (moldesign.Molecule): molecule to assign the forcefield to. It's ``ff``
               attribute will be overwritten with the new parameters
            display (bool): in jupyter, creates a parameterization display if the assignment fails
               (ignored outside of jupyter)

        Returns:
            Any: ``True`` if the forcefield can be assigned, ``ForcefieldAssignmentError``
               object otherwise
        """
        raise NotImplementedError()

    def add_ff(self, ff):
        """ Adds additional forcefield terms to this forcefield object.

        Note:
            These parameters are NOT currently checked for overlaps or compatibility with
            the pre-existing defintion!

        Args:
            ff (Forcefield): another forcefield to add
        """
        raise NotImplementedError()


@exports
class TLeapForcefield(Forcefield):
    """ Amber-type forcefield.

    Assignment routines for these forcefields are backed by ambertools (mostly tleap)

    Args:
        fflines (List[str]): commands to load this forcefield
        file_list (Dict[str, pyccc.FileReference]): any files necessary for assigning parameters
             (e.g. lib files, frcmod files)
    """
    TYPE = 'amber'

    def __str__(self):
        return "TLeap data at %x" % id(self)

    def __init__(self, fflines, file_list=None):
        self.names = None
        self._fflines = fflines
        self._file_list = file_list if file_list is not None else {}
        super().__init__()

    def assign(self, mol):
        newmol = self.create_prepped_molecule(mol)
        if newmol.num_atoms != mol.num_atoms:
            # TODO: much better error message here
            raise ForcefieldAssignmentError()

        else:
            # TODO: much more rigorous consistency checking
            newmol.ff.copy_to(mol)

    def create_prepped_molecule(self, mol, display=True):
        from ..interfaces import ambertools

        clean_molecule = ambertools._prep_for_tleap(mol)

        job = ambertools._run_tleap_assignment(clean_molecule, self._fflines, self._file_list)

        if 'output.inpcrd' in job.get_output():
            prmtop = job.get_output('output.prmtop')
            inpcrd = job.get_output('output.inpcrd')
            params = ambertools.AmberParameters(prmtop, inpcrd, job)
            m = mdt.read_amber(params.prmtop, params.inpcrd)
            newmol = mdt.helpers.restore_topology(m, mol)
            newmol.ff = mdt.forcefields.ForcefieldParams(newmol, params)
        else:
            newmol = None

        errors = ambertools._parse_tleap_errors(job, clean_molecule)

        show_parameterization_results(errors, clean_molecule, molout=newmol)

        if newmol is not None:
            return newmol
        else:
            raise ForcefieldAssignmentError(
                    'TLeap failed to assign force field parameters for %s' % mol, job)

    def add_ff(self, ff):
        self._fflines.extend(ff._fflines)
        for fname, fobj in ff._file_list.items():
            if fname in self._file_list:
                raise ValueError("Can't combine forcefields - two files with same name (%s)"
                                 %fname)
        self._file_list.update(ff._file_list)
