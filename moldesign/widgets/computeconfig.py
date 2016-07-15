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
import io
import os
import base64

import collections
import ipywidgets as ipy

import moldesign as mdt
from moldesign import compute, uibase

from . import toplevel
from . import __all__ as __packageall


@toplevel
def configure():
    from IPython.display import display
    display(MDTConfig())

# Some synonyms
about = configure
__packageall.append('about')


class MDTConfig(ipy.Box):
    def __init__(self):
        super(MDTConfig, self).__init__(display='flex', flex_flow='column')

        self.compute_config = ComputeConfig()
        self.tab_list = uibase.StyledTab([ipy.Box(),self.compute_config])
        self.tab_list.set_title(0, '^')
        self.tab_list.set_title(1, 'Compute configuration')

        self.children = [self.make_header(), self.tab_list]

    def make_header(self):
        img = io.open(os.path.join(mdt.PACKAGEPATH, '_static_data/img/banner.png'), 'r+b').read()
        encoded = base64.b64encode(img)
        img = '<img style="max-width:100%" src=data:image/png;base64,'+('%s>'%encoded)
        links = [self._makelink(*args) for args in
                   (("http://moldesign.bionano.autodesk.com/", 'About'),
                    ("https://forum.bionano.autodesk.com/c/Molecular-Design-Toolkit", 'Forum'),
                    ("https://github.com/autodesk/molecular-design-toolkit/issues", 'Issues'),
                    ("http://bionano.autodesk.com/MolecularDesignToolkit/explore.html",
                     "Tutorials"),
                    ('http://autodesk.github.io/molecular-design-toolkit/', 'Documentation'),
                    )]
        linkbar = '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'.join(links)
        return ipy.HTML(("<span style='float:left;font-size:0.8em;font-weight:bold'>Version: "
                         "{version}</span>"
                         "<span style='float:right'>{linkbar}</span>"
                         "<p>{img}</p>").format(img=img, linkbar=linkbar, version=mdt.__version__))

    @staticmethod
    def _makelink(url, text):
        return '<a href="{url}" target="_blank" title="{text}">{text}</a>'.format(url=url,
                                                                                  text=text)


class ComputeConfig(ipy.Box):
    def __init__(self):
        super(ComputeConfig, self).__init__(display='flex', flex_flow='column')

        self.engine_dropdown = ipy.Dropdown(description='Compute engine',
                                            options=ENGINE_DISPLAY,
                                            value=ENGINES.keys()[0],
                                            height='30px')
        self.engine_dropdown.observe(self.update_engine_display)

        self.engine_config_description = ipy.HTML('description')
        self.engine_config_value = ipy.Text('blank')
        self.engine_config = ipy.HBox([self.engine_config_description,
                                       self.engine_config_value])

        self._reset_config_button = ipy.Button(description='Reset',
                                               tooltip='Reset to current configuration')
        self._apply_changes_button = ipy.Button(description='Apply',
                                                tooltip='Apply for this session')
        self._save_changes_button = ipy.Button(description='Make default',
                                               tooltip='Make this the default for new sessions')
        self._reset_config_button.on_click(self.reset_config)
        self._apply_changes_button.on_click(self.apply_config)
        self._save_changes_button.on_click(self.save_config)

        self.children = [self.engine_dropdown,
                         ipy.HTML('<hr>'),
                         ipy.HBox([self.engine_config_description,
                                   self.engine_config_value]),
                         ipy.HBox([self._reset_config_button,
                                   self._apply_changes_button,
                                   self._save_changes_button])
                         ]
        self.reset_config()

    def update_engine_display(self, *args):
        self.engine_config_value.disabled = False
        enginename = self.engine_dropdown.value
        enginespec = ENGINES[enginename]

        self.engine_config_description.value = enginespec['hostdescription'] + ':'
        if enginename == 'free-compute-cannon':
            self.engine_config_value.value = compute.FREE_COMPUTE_CANNON
            self.engine_config_value.disabled = True
        else:
            self.engine_config_value.value = compute.config[enginespec['configkey']]

    def reset_config(self, *args):
        """ Reset configuration in UI widget to the stored values
        """
        if self.engine_dropdown.value != compute.config.engine_type:
            self.engine_dropdown.value = compute.config.engine_type
        else:
            self.update_engine_display()

    def apply_config(self, *args):
        enginename = self.engine_dropdown.value
        compute.config.engine_type = enginename
        compute.config[ENGINES[enginename]['configkey']] = self.engine_config_value.value
        compute.reset_from_config()

    def save_config(self, *args):
        raise NotImplementedError


ENGINES = collections.OrderedDict((
    ('free-compute-cannon', {'displayname': "Public CloudComputeCannon Demo",
                             'hostdescription': 'Autodesk-sponsored cloud compute server',
                             'configkey':'default_ccc_host'}),

    ('cloud-compute-cannon', {'displayname': 'CloudComputeCannon',
                              'hostdescription': 'Server name with port (e.g.,   '
                                                 '"192.168.0.1:8000")',
                              'configkey':'default_ccc_host'}),

    ('docker', {'displayname': 'Docker',
                'hostdescription': 'Docker host with port (e.g., "localhost:2375")',
                'configkey':'default_docker_host'}),

    ('docker-machine', {'displayname': 'Docker Machine',
                        'hostdescription': 'Name of docker-machine (e.g., "default")',
                        'configkey':'default_docker_machine'}),
))

ENGINE_DISPLAY = collections.OrderedDict((v['displayname'],k) for k,v in ENGINES.iteritems())