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
import base64
import io
import logging
import os
from collections import OrderedDict

import IPython.display
import ipywidgets as ipy
import traitlets

import moldesign as mdt

from moldesign.widgets.components import StyledTab

STANDARD = 25  # logging level between INFO and WARN

# We'll use standard + extra logging levels:
# debug, info, *status*, *handled*, warning, error

root = logging.getLogger('moldesign')
root.setLevel(STANDARD)

# TODO: we need to handle logging outside the widget context - what if user is in CLI?

_prev_tabs = None
_current_tabs = None
_capture_enabled = False


class NotebookHeader(ipy.Box):
    def __init__(self):
        super(NotebookHeader, self).__init__()
        self.make_header()
        self.header = self.make_header()
        self.configuration = self.make_configurator()
        self.children = [self.header, self.configuration]

    def make_header(self):
        img = io.open(os.path.join(mdt.PACKAGEPATH, 'notebooks/img/Top.png'), 'r+b').read()
        encoded = base64.b64encode(img)
        html = ['<img style="max-width:100%" src=data:image/png;base64,'+('%s>'%encoded)]
        linkbar = [self._makelink(*args) for args in
                   (("http://moldesign.bionano.autodesk.com/", 'About'),
                    ("http://bionano.autodesk.com/MolecularDesignToolkit/explore.html",
                     "Tutorials"),
                    ('http://moldesign.readthedocs.io/', 'Documentation'),
                    ("https://forum.bionano.autodesk.com/c/Molecular-Design-Toolkit", 'Forum'),
                    ("https://github.com/autodesk/molecular-design-toolkit/issues", 'Issues'))]
        html.append('&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'.join(linkbar))
        return ipy.HTML('\n'.join(html))

    def make_configurator(self):
        self._engine = ipy.Dropdown(description='Compute engine',
                                 options=['CloudComputeCannon',
                                          'LocalComputeCannon',
                                          'Docker',
                                          'Local process'])
        self._host = ipy.Text(description='Host')
        self._repo= ipy.Text(description='Docker image repository')
        self._version = ipy.Text(description='Docker version tag')
        self._restore_default_button = ipy.Button(description='Restore defaults')
        self._reset_config_button = ipy.Button(description='Reset', onclick=self.reset_config)
        self._apply_changes_button = ipy.Button(description='Apply')

        cfger = ipy.Box([self._host, self._repo, self._version,
                         self._restore_default_button, self._reset_config_button,
                         self._apply_changes_button],
                        layout=ipy.Layout(display='flex', flex_flow='row wrap'))
        return cfger

    def reset_config(self): pass


    @staticmethod
    def _makelink(url, text):
        return '<a href="{url}" target="_blank" title="{text}">{text}</a>'.format(url=url,
                                                                                  text=text)


def display(obj, title=None):
    """
    Registers a new view. This is mostly so that we can
    display all views from a cell in a LoggingTabs object.
    :param obj: The object to display. If it has a "get_display_object" method, \
        its return value is displayed
    :param title: A name for the object (otherwise, str(obj) is used)
    :return:
    """
    # for compatibility with earlier versions
    if hasattr(obj, 'get_display_object'):
        obj = obj.get_display_object()

    if _current_tabs is None:  # just display the damn thing
        IPython.display.display(obj)
        return

    _current_tabs.add_display(obj, title=title, display=True)


class WidgetValueHandler(logging.Handler):
    """
    Appends logging messages to a widget's value field.
    Watch this for performance - we're repeatedly +='ing a string, which is not good
    """
    def __init__(self, widget=None):
        """
        Initialize the handler.
        If stream is not specified, sys.stderr is used.
        """
        super(WidgetValueHandler, self).__init__()
        self.buffer = []
        if widget is None:
            self.widget = ipy.Textarea()
        else:
            self.widget = widget

    def emit(self, record):
        """
        Emit a record.
        If a formatter is specified, it is used to format the record.
        The record is then written to the stream with a trailing newline.  If
        exception information is present, it is formatted using
        traceback.print_exception and appended to the stream.  If the stream
        has an 'encoding' attribute, it is used to determine how to do the
        output to the stream.
        """
        try:
            msg = self.format(record)
            self.widget.value += msg
        except:
            self.handleError(record)


def enable_logging_widgets(enable=True):
    """

    :param enable: if False, disable it
    :return:
    """
    global _capture_enabled, _current_tabs

    from IPython import get_ipython
    ip = get_ipython()
    ":type: IPython.core.interactiveshell.InteractiveShell"

    if ip is None: return  # this isn't an interactive session

    if enable and not _capture_enabled:
        _capture_enabled = True
        #root.addHandler(capture_handler)
        ip.events.register('pre_run_cell', _capture_logging_displays)
        ip.events.register('post_run_cell', _finalize_logging_displays)

    elif (not enable) and _capture_enabled:
        _capture_enabled = False
        _current_tabs = None
        #root.removeHandler(capture_handler)
        ip.events.unregister('pre_run_cell', _capture_logging_displays)
        ip.events.unregister('post_run_cell', _finalize_logging_displays)

# TODO: Remove this workaround for https://github.com/ipython/ipywidgets/issues/628
# If ipywidgets#628 is a bug, then this can just be remoed
# If this is the expected behavior, we should come up with something more elegant
try:
    ipy.Text()
except traitlets.TraitError:
    widgets_enabled = False
else:
    widgets_enabled = True

if widgets_enabled:
    class LoggingTabs(StyledTab):
        def __init__(self, objects, display=False, **kwargs):
            """
            :param objects: dict of form {TITLE: <display object>}
            :param display: directly display the display collection
            :param kwargs: kwargs to pass to ipywidgets initializers
            """
            self.objs = OrderedDict(objects)
            super(LoggingTabs, self).__init__(objects.values(), **kwargs)
            self.selected_index = -1
            for ikey, key in enumerate(objects.iterkeys()):
                self.set_title(ikey, key)
            self._displayed = False
            if display:
                self._displayed = True
                IPython.display.display(self)

        def add_display(self, obj, title=None, display=True):
            title = mdt.utils.if_not_none(title, str(obj))
            title = title[:40]
            oldtitle = title
            ititle = 0
            while title in self.objs:
                ititle += 1
                title = '%s.%d' % (oldtitle, ititle)
            self.objs[title] = obj
            self.children += (obj,)  # ipywidgets requires a tuple for some reason
            self.set_title(len(self.children) - 1, title)
            if display and not self._displayed:
                IPython.display.display(self)
                self._displayed = True


class Logger(ipy.Textarea if widgets_enabled else object):
    # TODO: need a javascript-side widget that will accept a stream - string concatenation is bad
    def __init__(self, title='log', **kwargs):
        kwargs.setdefault('width', 400)
        kwargs.setdefault('height', 300)
        self.title = title
        if widgets_enabled:  # try to intialize widget
            super(Logger, self).__init__(**kwargs)
            self._is_widget = True
        else:
            self._is_widget = False
        self.active = False
        self.disabled = True  # so user can't overwrite

    def _write(self, string):
        if self._is_widget:
            if not self.active:
                display(self, self.title)
                self.active = True
            self.value += string.strip() + '\n'
        else:
            print string.strip()

    # temporary so that we can use this like a logging module later
    error = warning = info = handled = debug = status = _write


def _capture_logging_displays(display=False, **kwargs):
    global _current_tabs, _prev_tabs
    _prev_tabs = _current_tabs
    if widgets_enabled:
        _current_tabs = LoggingTabs(OrderedDict(x=ipy.Box()), display=display, **kwargs)
    else:
        _current_tabs = None
        enable_logging_widgets(False)
        print 'Failed to create UI logging system. Logging widgets disabled'


def _finalize_logging_displays(display=True, **kwargs):
    pass

# FOR NOW, *always* enable the logging widgets
enable_logging_widgets(True)
