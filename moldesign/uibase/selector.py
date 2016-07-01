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
# TODO: catch and log event exceptions
import ipywidgets as ipy

import moldesign as mdt
from moldesign import utils, viewer


class Selector(object):
    """
    This is the abstract base class for something that can make a selection.
    """

    def __init__(self, *args, **kwargs):
        self.selection_group = None
        self.selection_id = None
        super(Selector, self).__init__(*args, **kwargs)

    def handle_selection_event(self, selection):
        raise NotImplementedError()

    def fire_selection_event(self, new_selection):
        self.selection_group.update_selections(self, new_selection)


class SelectionGroup(ipy.Box):
    """
    Broadcasts selections among a group of widgets.
    It doesn't do much beside rebroadcast events:
    A SelectionGroup object will call the  "handle_selection_event" methods of
    all of its children whenever its "update_selections" method is called.
    """

    def __reduce__(self):
        """These don't gat passed around,
        so it reduces to NOTHING"""
        return utils.make_none, tuple()

    def __init__(self, *args, **kwargs):
        super(SelectionGroup, self).__init__(*args, **kwargs)
        self.selection_listeners = []
        self.num_listeners = 0
        self.selection = {}
        self.viewer = None
        self.graphviewer = None
        self.register_selection_listeners()

    def update_selections(self, event_source, selection):
        self.selection.update(selection)
        for listener in self.selection_listeners:
            listener.handle_selection_event(selection)

    def register_selection_listeners(self):
        self.num_listeners = 0
        self.selection_listeners = self.get_child_listeners(self)

    def get_child_listeners(self, element):
        if hasattr(element, 'handle_selection_event'):
            self.num_listeners += 1
            listeners = [element]
            element.selection_group = self
            element.selection_id = self.num_listeners
            if issubclass(element.__class__, mdt.viewer.GeometryViewer):
                self.viewer = element
            if issubclass(element.__class__, mdt.viewer.ChemicalGraphViewer):
                self.graphviewer = element
        else:
            listeners = []

        if hasattr(element, 'children'):
            for child in element.children:
                listeners.extend(self.get_child_listeners(child))
        return listeners

    @utils.args_from(viewer.GeometryViewer.set_color)
    def set_color(self, *args, **kwargs):
        if self.graphviewer: self.graphviewer.set_color(*args, **kwargs)
        if self.viewer: self.viewer.set_color(*args, **kwargs)

    @utils.args_from(viewer.GeometryViewer.set_color)
    def color_by(self, *args, **kwargs):
        if self.graphviewer: self.graphviewer.color_by(*args, **kwargs)
        if self.viewer: self.viewer.color_by(*args, **kwargs)

    @utils.args_from(viewer.GeometryViewer.set_color)
    def set_colors(self, *args, **kwargs):
        if self.graphviewer: self.graphviewer.set_colors(*args, **kwargs)
        if self.viewer: self.viewer.set_colors(*args, **kwargs)

    @utils.args_from(viewer.GeometryViewer.unset_color)
    def unset_color(self, *args, **kwargs):
        if self.graphviewer: self.graphviewer.unset_color(*args, **kwargs)
        if self.viewer: self.viewer.unset_color(*args, **kwargs)

    def __getattr__(self, item):
        if self.viewer is not None: return getattr(self.viewer, item)
        else: raise AttributeError(item)


class ValueSelector(Selector):
    """ This is an abstract mixin for a widget class
    """
    def __init__(self, value_selects=None, **kwargs):
        self.value_selects = value_selects
        super(ValueSelector, self).__init__(**kwargs)
        self.observe(self.value_update, 'value')
        self.__hold_fire = False

    def value_update(self, *args):
        if self.__hold_fire: return  # prevent recursive selections
        self.fire_selection_event({self.value_selects: self.value})

    def handle_selection_event(self, selection):
        self.__hold_fire = True
        try:
            if self.value_selects in selection:
                self.value = selection[self.value_selects]
        except Exception as exc:
            print 'ERROR: (ignored) %s' % exc
        self.__hold_fire = False


def create_value_selector(widget, value_selects, **kwargs):
    """
    Creates a UI element (slider, checkbox, etc.) to add to your
    selection group.
    :param widget: widget class (e.g. ipywidgets.FloatSlider)
    :param value_selects: What the value of the widget selects (e.g. time)
    :param kwargs: keyword arguments for the widget (e.g. description = "time")
    :return: the constructed widget
    """

    class SelectionWidget(ValueSelector, widget): pass

    return SelectionWidget(value_selects=value_selects, **kwargs)


