import os
import yaml
import pytest
import moldesign as mdt
from moldesign.compute import packages

from .molecule_fixtures import *

__PYTEST_MARK__ = ['internal', 'config']
# mark all tests in this module with this label (see ./conftest.py)


@pytest.fixture(scope='function')
def tempconfigfile(tmpdir):
    tmpdir = str(tmpdir)
    path = os.path.join(tmpdir, 'moldesign.yml')
    oldenviron = os.environ.get('MOLDESIGN_CONFIG', None)
    os.environ['MOLDESIGN_CONFIG'] = path
    with open(path, 'w') as cfile:
        cfile.write('devmode: true')

    yield path

    if oldenviron is not None:
        os.environ['MOLDESIGN_CONFIG'] = oldenviron
    else:
        os.environ.pop('MOLDESIGN_CONFIG')
    clean_config()


def clean_config():
    mdt.compute.config.clear()
    mdt.compute.init_config()


def test_write_conf(tempconfigfile):
    mdt.compute.update_saved_config(run_local={'nwchem':True, 'opsin':True})
    mdt.compute.update_saved_config(run_remote={'openmm':True}, run_local={'nwchem':False})
    mdt.compute.update_saved_config(devmode=True)

    with open(tempconfigfile, 'r') as conffile:
        data = yaml.load(conffile)

    assert data['run_local'] == {'nwchem':False, 'opsin':True}
    assert data['run_remote'] == {'openmm':True}
    assert data['devmode']


@pytest.fixture
def disconnect_docker():
    oldhost = mdt.compute.config.default_docker_host
    mdt.compute.config.default_docker_host = '000.111:333'
    yield
    mdt.compute.config.default_docker_host = oldhost
    mdt.compute.compute.default_engine = None


def test_docker_connection_failure(disconnect_docker):
    with pytest.raises(mdt.exceptions.DockerError):
        mdt.compute.reset_compute_engine()


@pytest.mark.skipif(not packages.openbabel.is_installed(),
                    reason='Only run if openbabel is present locally')
def test_set_docker_python(tempconfigfile):
    njob = mdt._njobs

    g1 = mdt.from_smiles('CC')

    assert mdt._njobs == njob

    mdt.compute.update_saved_config(run_remote={'openbabel':True})
    clean_config()

    g2 = mdt.from_smiles('CC')
    assert mdt._njobs == njob + 1


def test_docker_image_strings(tempconfigfile):
    mdt.compute.config.devmode = True
    assert mdt.compute.get_image_path('blah') == 'blah:dev'
    assert mdt.compute.get_image_path('a/b') == 'a/b:dev'

    mdt.compute.config.devmode = False
    mdt.compute.config.default_repository = 'myrepo/thing:'
    mdt.compute.config.default_version_tag = 'tagtag'
    assert mdt.compute.get_image_path('blah') == 'myrepo/thing:blah-tagtag'

    mdt.compute.config.default_repository = 'myrepo/thing'
    mdt.compute.config.default_version_tag = 'tagtag'
    assert mdt.compute.get_image_path('blah') == 'myrepo/thing/blah:tagtag'

    mdt.compute.config.default_repository = 'myrepo/thing/'
    assert mdt.compute.get_image_path('blah') == 'myrepo/thing/blah:tagtag'


def test_nbmolviz_errors(ethylene):
    """ These should always raise importerrors because we're not running in Jupyter
    """
    with pytest.raises(ImportError):
        mdt.widgets.BondSelector()

    with pytest.raises(ImportError):
        ethylene.draw()
