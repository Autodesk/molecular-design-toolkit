import collections

import ipywidgets as ipy

from moldesign.uibase import UnitText


class Configurator(ipy.Box):
    """ Interactive configuration widget for a list of parameters

    Args:
        paramlist (Mapping[str,object]): Mapping of parameter names to their current values
        paramdefs (Mapping[str,Parameter]): Mapping of parameter names to their definitions

    Attributes:
        selectors (Mapping[str,ParamSelector]): Mapping of parameters names to their selector
           widgets
    """

    def __init__(self, paramlist, paramdefs):
        super(Configurator, self).__init__(layout=ipy.Layout(display='flex', flex_flow='column'))
        self.paramlist = paramlist

        self.selectors = collections.OrderedDict([(p.name, ParamSelector(p)) for p in paramdefs])
        self.reset_values()

        self.children = self.selectors

    def reset_values(self):
        reset_params = set()
        for name, value in self.paramlist.iteritems():
            self.selectors[name].selector.value = value
            reset_params.add(name)

        for name, paramdef in self.paramdefs.iteritems():
            if name not in reset_params:  # if it's not set already, make it the default
                self.selectors[name].default()

    def apply_values(self):
        for paramname, selector in self.selectors.iteritems():
            self.paramlist[paramname] = selector.selector.value


class ParamSelector(ipy.Box):
    def __init__(self, param):
        super(ParamSelector, self).__init__(layout=ipy.Layout(display='flex', flex_flow='row wrap'))

        self.param = param

        children = []
        self.name = ipy.HTML(param.name)
        children.append(self.name)

        if param.choices:
            self.selector = ipy.Dropdown(options=param.choices)
        elif param.type == bool:
            self.selector = ipy.ToggleButtons(options=[True, False])
        elif param.units:
            self.selector = UnitText(units=param.units)
        elif param.type == float:
            self.selector = ipy.FloatText()
        elif param.type == int:
            self.selector = ipy.IntText()
        elif param.type == str:
            self.selector = ipy.Text()
        else:
            self.selector = ipy.HTML('<i>unknown datatype</i>')
        children.append(self.selector)

        children = [self.name, self.selector]

        self.default_button = None
        if param.default:
            self.default_button = ipy.Button(description='Default',
                                             tooltip='Reset to: %s' % self.param.default,
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