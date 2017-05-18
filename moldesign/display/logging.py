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


try:
    from nbmolviz.uibase import Logger, display_log
except ImportError:
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
            print string.strip()

        # temporary so that we can use this like a logging module later
        error = warning = info = handled = debug = status = _write

    def display_log(obj, title=None, show=False):
        """
        Registers a new view. This is mostly so that we can
        display all views from a cell in a LoggingTabs object.
        :param obj: The object to display. If it has a "get_display_object" method, \
            its return value is displayed
        :param title: A name for the object (otherwise, str(obj) is used)
        :return:
        """
        print obj