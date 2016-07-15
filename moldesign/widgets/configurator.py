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

import collections

import ipywidgets as ipy

from moldesign.uibase import UnitText
from moldesign import utils


def exports(o):
    __all__.append(o.__name__)
    return o
__all__ = []


@exports
class Configurator(ipy.Box):
    """ Interactive configuration widget for a list of parameters

    Args:
        paramlist (Mapping[str,object]): Mapping of parameter names to their current values
        paramdefs (Mapping[str,Parameter]): Mapping of parameter names to their definitions

    Attributes:
        selectors (Mapping[str,ParamSelector]): Mapping of parameters names to their selector
           widgets
    """
    # TODO: help text. Use javascript to show/hide? http://stackoverflow.com/a/1313353/1958900

    def __init__(self, paramlist, paramdefs, title=None):
        super(Configurator, self).__init__(layout=ipy.Layout(display='flex',
                                                             flex_flow='column',
                                                             align_self='flex-start',
                                                             align_items='stretch',
                                                             max_width='100%'))
        self.paramlist = paramlist
        self.paramdefs = paramdefs

        self.apply_button = ipy.Button(description='Apply')
        self.apply_button.on_click(self.apply_values)

        self.reset_button = ipy.Button(description='Reset')
        self.reset_button.on_click(self.reset_values)
        self.buttons = ipy.Box([self.reset_button, self.apply_button],
                               layout=ipy.Layout(align_self='center'))

        self.selectors = collections.OrderedDict([(p.name, ParamSelector(p)) for p in paramdefs])
        self.reset_values()

        title = utils.if_not_none(title, 'Configuration')
        self.title = ipy.HTML('<center><h4>%s</h4></center><hr>' % title,
                              align_self='center')

        self.currentconfig = ipy.Textarea(description='Current params:',
                                          disabled=True,
                                          value=str(paramlist).replace(', ', ',\n   '),
                                          width='350px')
        self.middle = ipy.HBox([ipy.VBox(self.selectors.values()), self.currentconfig])
        self.children = [self.title, self.middle, self.buttons]

    def reset_values(self, *args):
        reset_params = set()
        for name, value in self.paramlist.iteritems():
            self.selectors[name].selector.value = value
            reset_params.add(name)

        for paramdef in self.paramdefs:
            if paramdef.name not in reset_params:  # if it's not set already, make it the default
                self.selectors[paramdef.name].default()

    def apply_values(self, *args):
        for paramname, selector in self.selectors.iteritems():
            self.paramlist[paramname] = selector.selector.value
        self.currentconfig.value = str(self.paramlist).replace(', ', ',\n   ')


class ParamSelector(ipy.Box):

    WIDGETKWARGS = {'width':'200px'}

    def __init__(self, param):
        super(ParamSelector, self).__init__(layout=ipy.Layout(display='flex',
                                                              flex_flow='nowrap',
                                                              align_content='baseline'))

        self.param = param

        children = []
        self.name = ipy.HTML("<p style='text-align:right'>%s:</p>" % param.displayname,
                             width='200px')
        children.append(self.name)

        if param.choices:
            self.selector = ipy.Dropdown(options=param.choices, **self.WIDGETKWARGS)
        elif param.type == bool:
            self.selector = ipy.ToggleButtons(options=[True, False], **self.WIDGETKWARGS)
        elif param.units:
            self.selector = UnitText(units=param.units, **self.WIDGETKWARGS)
        elif param.type == float:
            self.selector = ipy.FloatText(**self.WIDGETKWARGS)
        elif param.type == int:
            self.selector = ipy.IntText(**self.WIDGETKWARGS)
        elif param.type == str:
            self.selector = ipy.Text(**self.WIDGETKWARGS)
        else:
            self.selector = ipy.HTML('<i>unknown datatype</i>')
        children.append(self.selector)

        children = [self.name, self.selector]

        self.default_button = None
        if param.default:
            self.default_button = ipy.Button(description='Default',
                                             tooltip='Set to default: %s' % self.param.default,
                                             width='75px')
            self.default_button.on_click(self.default)
            children.append(self.default_button)
            self.default()

        self.help_link = None
        if param.help_url:
            self.help_link = ipy.HTML('<a href="%s" target="_blank">?</a>' % param.help_url)
            children.append(self.help_link)

        self.children = children

    def default(self, *args):
        self.selector.value = self.param.default

    @property
    def value(self):
        return self.selector.value

    @value.setter
    def value(self, v):
        self.selector.value = v