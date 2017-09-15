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

from ..utils import exports_names, exports
from ..widgets import nbmolviz_enabled
from .. import units as u


if nbmolviz_enabled:
    from nbmolviz.uielements import Logger, display_log
    exports_names('Logger', 'display_log')

else:
    @exports
    class Logger(object):
        def __init__(self, title='log', **kwargs):
            kwargs.setdefault('width', '75%')
            kwargs.setdefault('height', '300px')
            kwargs.setdefault('font_family', 'monospace')
            self.title = title
            self._is_widget = False
            self.active = False
            self.disabled = True  # so user can't overwrite

        def _write(self, string):
            print(string.strip())

        # temporary so that we can use this like a logging module later
        error = warning = info = handled = debug = status = _write

    @exports
    def display_log(obj, title=None, show=False):
        """
        Registers a new view. This is mostly so that we can
        display all views from a cell in a LoggingTabs object.
        :param obj: The object to display. If it has a "get_display_object" method, \
            its return value is displayed
        :param title: A name for the object (otherwise, str(obj) is used)
        :return:
        """
        print(obj)


@exports
class DynamicsLog(object):
    ROW_FORMAT = ("{:<10.2f}") + 3*(" {:>15.4f}")
    HEADER_FORMAT = ROW_FORMAT.replace('.4f','s').replace('.2f','s')

    def __init__(self):
        self._printed_header = False

    def print_header(self):
        timeheader = 'time /'
        peheader = 'potential /'
        keheader = 'kinetic /'
        temperatureheader = 'T /'
        print(self.HEADER_FORMAT.format(timeheader, peheader, keheader, temperatureheader))

        timeunits = '{}'.format(u.default.time)
        peunits = '{}'.format(u.default.energy)
        keunits = '{}'.format(u.default.energy)
        temperatureunits = '{}'.format(u.default.temperature)
        print(self.HEADER_FORMAT.format(timeunits, peunits, keunits, temperatureunits))

        self._printed_header = True

    def print_step(self, mol, properties):
        from . import kinetic_energy, kinetic_temperature
        if not self._printed_header:
            self.print_header()

        ke = kinetic_energy(properties['momenta'], mol.masses)
        t = kinetic_temperature(ke, mol.dynamic_dof)
        print(self.ROW_FORMAT.format(properties['time'].defunits_value(),
                                     properties['potential_energy'].defunits_value(),
                                     ke.defunits_value(),
                                     t.defunits_value()))
