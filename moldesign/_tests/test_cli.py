import subprocess
import pytest

import moldesign as mdt
from moldesign.external import pathlib


@pytest.fixture
def example_path(tmpdir):
    path = pathlib.Path(str(tmpdir))
    subprocess.check_call('python -m moldesign copyexamples',
                          shell=True,
                          cwd=str(path))
    return path


def test_exampled_copied(example_path):
    path = example_path
    assert (path / 'moldesign-examples').is_dir()
    with (path / 'moldesign-examples' / '.mdtversion').open('r') as verfile:
        assert verfile.read().strip() == mdt.__version__


def test_no_overwrite_examples(example_path):
    path = example_path
    try:
        subprocess.check_call('python -m moldesign copyexamples',
                              shell=True,
                              cwd=str(path))
    except subprocess.CalledProcessError as err:
        assert err.returncode == 200

    else:
        assert False, "Expected CalledProcessError"


def test_example_version_warning(tmpdir):
    path = example_path(tmpdir)  # call this directly because we mangle the test dir version
    with (path / 'moldesign-examples' / '.mdtversion').open('w') as verfile:
        verfile.write(u'0.1.0')

    try:
        subprocess.check_call('python -m moldesign copyexamples',
                              shell=True,
                              cwd=str(path))
    except subprocess.CalledProcessError as err:
        assert err.returncode == 201

    else:
        assert False, "Expected CalledProcessError"


def test_version_command():
    ver = subprocess.check_output('python -m moldesign version'.split()).splitlines()[-1]
    assert ver.decode('ascii') == mdt.__version__


def test_dumpenv_command():
    # just test that it doesn't error
    subprocess.check_call('python -m moldesign dumpenv'.split())


def test_print_environment():
    # just test that it still works
    mdt.data.print_environment()

def test_config_command():
    # just test that it doesn't error
    subprocess.check_call('python -m moldesign config'.split())
