from __future__ import print_function, absolute_import, division
from future.builtins import *
from future import standard_library
standard_library.install_aliases()

# Copyright 2017 Autodesk Inc.
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
from .. import utils
from ..helpers import display_log


def get_image_path(image_name, _devmode=None):
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

    if _devmode is None:
        _devmode = config.devmode

    if _devmode:
        return image_name + ':dev'

    if not config.default_repository:
        name = image_name
    elif config.default_repository[-1] in '/:':
        name = '%s%s' % (config.default_repository, image_name)
    else:
        name = '%s/%s' % (config.default_repository, image_name)

    if not config.default_repository:
        img = name
    elif config.default_repository[-1] == ':':
        img = '%s-%s' % (name, config.default_version_tag)
    elif config.default_version_tag:
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



def run_job(job, engine=None, wait=True, jobname=None, display=True,
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
        display_log(job.get_display_object(), jobname)

    if wait:
        job.wait()
        if _return_result: return job.result

    return job


@utils.args_from(run_job, only='engine wait jobname display'.split())
def runremotely(func, args=None, kwargs=None,
                jobname=None, engine=None, image=None, wait=True, display=True,
                persist_refs=True, when_finished=None):
    """ Runs a python command remotely.

    Args:
        job (pyccc.Job): The job to run

    Returns:
        pyccc.PythonJob OR object: reference to the job if wait=False, or the function's
           return value, if wait=True
    """
    import pyccc

    if args is None:
        args = []
    if kwargs is None:
        kwargs = {}

    if image is None:
        image = mdt.compute.config.default_python_image

    if jobname is None:
        jobname = func.__name__
        if args:
            jobname += str(args[0])

    call = pyccc.PythonCall(func, *args, **kwargs)
    job = pyccc.PythonJob(command=call, image=image, engine=engine, name=jobname,
                          submit=False, persist_references=persist_refs,
                          sendsource=False, when_finished=when_finished)
    job = run_job(job, wait=wait, display=display)
    if wait:
        return job.result
    else:
        return job


