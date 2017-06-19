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
from .helpers.widgets import nbmolviz_installed

if nbmolviz_installed:
    from nbmolviz.widgets import BondSelector, GeometryBuilder, ResidueSelector, Symmetrizer
    from nbmolviz.widgets.computeconfig import configure, about
else:
    from .helpers.widgets import not_installed_method
    BondSelector = GeometryBuilder = ResidueSelector = Symmetrizer = configure = about = \
        not_installed_method

__all__ = 'BondSelector GeometryBuilder ResidueSelector Symmetrizer'.split()
