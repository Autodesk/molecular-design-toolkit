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
import moldesign as mdt
from moldesign import utils


def get_image_path(image_name):
    """ Returns a fully qualified tag that points to the correct registry

    Args:
        image_name (str): name of the image (without tags, repository, etc.)

    Examples:
        >>> config.update({'default_repository':'my.docker.registry/orgname/myrepo:',
                           'default_version_tag':'latest'})
        >>> get_image_path('myimage')
        'my.docker.registry/orgname/myrepo:myimage-latest'

        >>> config.update({'default_repository':'docker.io/myorg',
                           'default_version_tag':'0.2'})
        >>> get_image_path('someimage')
        'docker.io/myorg/someimage:0.2
    """
    from . import config

    if not config.default_repository:
        name = image_name
    else:
        name = '%s%s' % (config.default_repository, image_name)

    if not config.default_repository:
        img = name
    elif config.default_repository[-1] == ':':
        img = '%s-%s' % (name, config.default_version_tag)
    elif config.version_tag:
        img = '%s:%s' % (name, config.default_version_tag)
    else:
        raise ValueError('Faulty docker repository configuration not recognized')

    return img


class DummyJob(object):
    """
    A job that doesn't actually need to run.
    Useful as a return value for processes that return a job object.
    """
    status = 'finished'
    wait = kill = lambda self: None
    get_stdout_stream = get_stderr_stream = lambda self: ''
    stdout = stderr = ''

    @staticmethod
    def get_output_files(filename=None):
        if filename is not None:
            return {}[filename]  # raises the exception on purpose
        else:
            return {}

    def __init__(self, result, updated_object=None):
        self.result = result
        if updated_object:
            self.updated_object = updated_object


def run_job(job, engine=None, image=None, wait=True, jobname=None, display=True,
            _return_result=False):
    """ Helper for running jobs.

    Args:
        job (pyccc.Job): The job to run
        engine (pyccc.Engine): Engine to run this job on (default:
            ``moldesign.compute.get_engine()``)
        image (str): URL for the docker image
        wait (bool): if True, block until this function completes and return the function's
            return value. Otherwise, return a job object immediately that can be queried later.
        display (bool): if True, show logging output for this job

    Returns:
        pyccc job object OR function's return value
    """

    engine = utils.if_not_none(engine, mdt.compute.get_engine())

    if engine is None:
        raise ValueError('No compute engine configured! Configure MDT using '
                         'moldesign.compute.config')

    engine.submit(job)

    jobname = utils.if_not_none(jobname, job.name)

    if display:
        mdt.uibase.display_log(job.get_display_object(), jobname)

    if wait:
        job.wait()
        if _return_result: return job.result

    return job






