"""
This module contains abstract base classes for potential models, integrators, and various
associated data types (force fields, orbitals, basis sets, etc.).
"""

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
from future.utils import with_metaclass
import funcsigs

import moldesign as mdt
from .helpers import WidgetMethod
from . import utils


class _InitKeywordMeta(type):
    """ Constructs a custom call signature for __init__ based on cls.PARAMETERS.
    """
    @property
    def __signature__(self):
        if hasattr(self, '__customsig'):
            return self.__customsig

        kwargs = []
        for param in self.PARAMETERS:
            kwargs.append(funcsigs.Parameter(param.name,
                                             default=param.default,
                                             kind=funcsigs.Parameter.POSITIONAL_OR_KEYWORD))

        self.__customsig = funcsigs.Signature(kwargs, __validate_parameters__=True)
        return self.__customsig


class Method(with_metaclass(_InitKeywordMeta, object)):
    """Abstract Base class for energy models, integrators, and "heavy duty" simulation objects

    Args:
        **kwargs (dict): list of parameters for the method.

    Attributes:
       mol (mdt.Molecule): the molecule this method is associated with
    """

    PARAMETERS = []
    """ list: list of Parameters that can be used to configure this method
    """

    PARAM_SUPPORT = {}
    """ Mapping(str, list): List of supported values for parameters (if a parameter is not found,
     it's assumed that all possible values are supported)
    """

    configure = WidgetMethod('method.configure')

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
        self.params = utils.DotDict(params)
        # Set default parameter values
        for param in self.PARAMETERS:
            if param.name not in self.params:
                self.params[param.name] = param.default

    def __eq__(self, other):
        return self.__class__ is other.__class__ and self.params == other.params

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