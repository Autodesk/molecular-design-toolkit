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
"""
This module contains abstract base classes for potential models, integrators, and various
associated data types (force fields, orbitals, basis sets, etc.).
"""

from moldesign.utils import DotDict


class Method(object):
    """Abstract Base class for energy models and integrators

    Args:
        **kwargs (dict): list of parameters for the method.

    Attributes:
       mol (mdt.Molecule): the molecule this method is associated with
    """

    PARAMETERS = []
    """ list: list of Parameters that can be used to configure this method
    """

    def __reduce__(self):
        return _make_method, (self.__class__, self.params, self.mol)

    def __init__(self, **params):
        """
        :param params:
        :return:
        """
        # TODO: better documentation for the expected keywords

        self._prepped = False
        self.status = None
        self.mol = None
        self.params = DotDict(params)
        # Set default parameter values
        for param in self.PARAMETERS:
            if param.name not in self.params:
                self.params[param.name] = param.default

    def configure(self):
        from moldesign.widgets.configurator import Configurator
        return Configurator(self.params, self.PARAMETERS,
                            title='Configuration for %s' % self.__class__.__name__)

    @classmethod
    def get_parameters(cls):
        """
        This doesn't do anything right now except provide guidelines for programmers
        """
        return cls.PARAMETERS

    def get_forcefield(self):
        raise NotImplementedError()

    @classmethod
    def print_parameters(cls):
        params = cls.PARAMETERS
        lines = []
        for obj in params:
            description = ''
            if obj.choices:
                description = '%s' % obj.choices
                if obj.types:
                    description += ' or '
            if obj.types:
                description += 'Type %s' % obj.types

            doc = '%s: %s (DEFAULT: %s)' % (obj.name, description, obj.default)
            lines.append(doc)
        return '\n'.join(lines)


def _make_method(cls, params, mol):
    """
    Helper for serialization - allows __reduce__ to use kwargs
    """
    obj = cls(**params)
    obj.mol = mol
    return obj