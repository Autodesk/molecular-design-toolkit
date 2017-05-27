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


class ConvergenceFailure(Exception):
    """ Raised when an iterative calculation fails to converge """
    pass


class NotCalculatedError(Exception):
    """ Raised when a molecular property is requested that hasn't been calculated """
    pass


class UnhandledValenceError(Exception):
    def __init__(self, atom):
        self.message = 'Atom %s has unhandled valence: %d' % (atom, atom.valence)


class QMConvergenceError(Exception):
    """ Raised when an iterative QM calculation (typically SCF) fails to converge
    """
    pass
