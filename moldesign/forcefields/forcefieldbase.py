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

from .exceptions import ForcefieldAssignmentError


class Forcefield(object):
    """ Abstract base class representing a biomolecular forcefield definition,
    such as amber14sb or charm22
    """

    TYPE = 'abstract'
    "str: FF type (e.g., amber, charmm, opls, polarizable)"

    FFNAME = 'abstract'
    "str: name of this forcefield (charm22, amber14sb, etc.)"


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

    def create_fixed_molecule(self, mol, display=True):
        """ Return a modify _copy_ of the passed molecule that's ready for forcefield assignment.

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


