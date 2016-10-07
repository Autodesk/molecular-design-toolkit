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
from pip._vendor.packaging import version

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
        self.changelog = ChangeLog()
        self.tab_list = uibase.StyledTab([ipy.Box(), self.compute_config, self.changelog])
        self.tab_list.set_title(0, '^')
        self.tab_list.set_title(1, 'Compute configuration')
        self.tab_list.set_title(2, "What's new")

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


class ChangeLog(ipy.Box):
    def __init__(self):
        super(ChangeLog, self).__init__(orientation='vertical')
        try:
            current = version.parse(mdt.__version__)
            latest = self.version_check()
            if current >= latest:
                versiontext = 'Up to date. Latest release: %s' % latest
            else:
                versiontext = ('New release available! '
                               '(Current: %s, latest: %s <br>' % (current, latest) +
                               '<b>Install it:</b> '
                               '<span style="font-family:monospace">pip install -U moldesign'
                               '</span>')
        except Exception as e:
            versiontext = '<b>Failed update check</b>: %s' % e

        self.version = ipy.HTML(versiontext)
        self.textarea = ipy.Textarea(width='700px', height='300px')

        p1 = os.path.join(mdt.PACKAGEPATH, "HISTORY.rst")
        p2 = os.path.join(mdt.PACKAGEPATH, "..", "HISTORY.rst")
        if os.path.exists(p1):
            path = p1
        elif os.path.exists(p2):
            path = p2
        else:
            path = None

        if path is not None:
            with open(path, 'r') as infile:
                self.textarea.value = infile.read()
        else:
            self.textarea.value = 'HISTORY.rst not found'

        self.textarea.disabled = True
        self.children = (self.version, self.textarea)

    @staticmethod
    def version_check():
        """
        References:
            http://code.activestate.com/recipes/577708-check-for-package-updates-on-pypi-works-best-in-pi/
        """
        import xmlrpclib
        pypi = xmlrpclib.ServerProxy('https://pypi.python.org/pypi')
        return pypi.package_releases('moldesign')


class ComputeConfig(ipy.Box):
    def __init__(self):
        super(ComputeConfig, self).__init__(display='flex', flex_flow='column')

        self.engine_dropdown = ipy.Dropdown(description='Compute engine',
                                            options=ENGINE_DISPLAY,
                                            value=ENGINES.keys()[0],
                                            height='30px')
        self.engine_dropdown.observe(self.update_engine_display)

        self.engine_config_description = ipy.HTML('description')
        self.engine_config_value = ipy.Text('blank', width='500px')
        self.engine_config = ipy.HBox([self.engine_config_description,
                                       self.engine_config_value])

        self._reset_config_button = ipy.Button(description='Reset',
                                               tooltip='Reset to current configuration')
        self._apply_changes_button = ipy.Button(description='Apply',
                                                tooltip='Apply for this session')
        self._save_changes_button = ipy.Button(description='Make default',
                                               tooltip='Make this the default for new sessions')
        self._test_button = ipy.Button(description='Test connection',
                                       tooltip='Test if MDT can run jobs with the '
                                               'current configuration')
        self._reset_config_button.on_click(self.reset_config)
        self._apply_changes_button.on_click(self.apply_config)
        self._save_changes_button.on_click(self.save_config)
        self._test_button.on_click(self.test_connection)

        self.children = [self.engine_dropdown,
                         ipy.HTML('<hr>'),
                         ipy.HBox([self.engine_config_description,
                                   self.engine_config_value]),
                         ipy.HBox([self._reset_config_button,
                                   self._apply_changes_button,
                                   self._test_button,
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
        if compute.config.engine_type not in ENGINES[self.engine_dropdown.value]['aliases']:
            self.engine_dropdown.value = compute.config.engine_type
        else:
            self.update_engine_display()

    def apply_config(self, *args):
        enginename = self.engine_dropdown.value
        compute.config.engine_type = ENGINES[enginename]['aliases'][0]
        compute.config[ENGINES[enginename]['configkey']] = self.engine_config_value.value
        compute.reset_compute_engine()

    def test_connection(self, *args):
        self.apply_config()
        engine = compute.default_engine
        if engine is None:
            raise ValueError('Failed to create compute engine with current configuration')
        engine.test_connection()
        print "SUCCESS: %s is accepting jobs" % engine

    def save_config(self, *args):
        compute.write_config()


class RegistryConfig(ipy.Box):
    def __init__(self):
        super(RegistryConfig, self).__init__(display='flex', flex_flow='column')
        self.repo_field = ipy.Text(description='Image repository')
        self.version_field = ipy.Text(description='Image version')

        self._reset_config_button = ipy.Button(description='Reset',
                                               tooltip='Reset to current configuration')

        self._apply_changes_button = ipy.Button(description='Apply',
                                                tooltip='Apply for this session')
        self._save_changes_button = ipy.Button(description='Make default',
                                               tooltip='Make this the default for new sessions')
        self._pull_button = ipy.Button(description='Pull images',
                                       tooltip=
                                       'Download all moldesign images to the compute engine')

        self.children = (ipy.HBox([self.repo_field, self.version_field]),
                         ipy.HBox([self._reset_config_button,
                                                self._apply_changes_button,
                                                self._pull_button]))

        self._reset_config_button.on_click(self.reset_config)
        self._apply_changes_button.on_click(self.apply_config)
        self._save_changes_button.on_click(self.save_config)
        self._test_button.on_click(self.test_connection)

    def reset_config(self, *args):
        self.repo_field.value = mdt.compute.config.default_repository
        self.version_field.value = mdt.compute.config.version_tag

    def apply_config(self, *args):
        compute.config.default_repository = self.repo_field.value
        compute.config.version_tag = self.version_field.value


_enginedefs = (
    ('free-compute-cannon', {'displayname': "Public CloudComputeCannon Demo",
                             'hostdescription': 'Autodesk-sponsored cloud compute server',
                             'configkey': 'default_ccc_host',
                             'aliases': ('ccc', 'cloudcomputecannon')
                             }),

    ('cloud-compute-cannon', {'displayname': 'CloudComputeCannon',
                              'hostdescription': 'Server address and port (e.g.,   '
                                                 '"192.168.0.1:9000")',
                              'configkey': 'default_ccc_host',
                              'aliases': ('ccc', 'cloudcomputecannon')}),

    ('docker', {'displayname': 'Docker',
                'hostdescription': 'Docker host with port (e.g., "localhost:2375")',
                'configkey': 'default_docker_host',
                'aliases': ('docker',)
                }),

    ('docker-machine', {'displayname': 'Docker Machine',
                        'hostdescription': 'Name of docker-machine (e.g., "default")',
                        'configkey': 'default_docker_machine',
                        'aliases': ('docker-machine',)
                        }),
)

ENGINES = collections.OrderedDict(_enginedefs)
ENGINE_DISPLAY = collections.OrderedDict((v['displayname'],k) for k,v in ENGINES.iteritems())