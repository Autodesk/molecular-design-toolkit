import pytest
import moldesign as mdt


def test_docker_image_strings():
    savedconfig = mdt.compute.config.copy()

    try:
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

    finally:
        mdt.compute.config = savedconfig



