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
import functools
import logging
import types

import moldesign as mdt
import pyccc.python as bpy
from moldesign.utils import DotDict, if_not_none
from moldesign import config, utils


class BuckyLog(logging.Logger):
    pass


def get_image_path(image_name, engine=None):
    """
    Returns a fully qualified tag that points to the correct registry
    (eventually)
    """
    if not config.cfg.default_repository:
        name = image_name
    else:
        name = '%s%s' % config.cfg.default_repository, image_name

    if not config.cfg.default_repository:
        img = name
    elif config.cfg.default_repository[-1] == ':':
        img = '%s-%s' % (name, config.cfg.version_tag)
    elif config.cfg.version_tag:
        img = '%s:%s' % (name, config.cfg.version_tag)

    # TODO: what happens if the engine also offers a "get_image_tag" (change that name) method?
    return img


launch = config.cfg.default_engine.launch  # this won't get updated, but it keeps the call signature
default_engine = config.cfg.default_engine
default_image = config.cfg.default_image


class RunsRemotely(object):
    def __init__(self, remote=True, display=True,
                 jobname=None,
                 sendsource=False,
                 engine=None,
                 image=None,
                 outputs=None,
                 cfg=None,
                 is_imethod=False):
        """Function decorator to run a job remotely.

        Note:
         This ONLY works for pure functions - where you're interested in the
           return value only. Side effects won't be visible to the user.

        Args:
            remote (bool): Whether or not to run the job remotely (for toggling remote functionality)
            display (bool): Create a jupyter logging display for the remote job
            jobname (str): Name metadata - defaults to the __name__ of the function
            sendsource (bool): if False (default), call this function directly on the remote worker;
               if True, send the function's source code (for debugging, mostly)
            engine (pyccc.engine.EngineBase): engine to send the job to
            image (str): name of the docker image (including registry, repository, and tags)
            is_imethod (bool): This is an instancemethod (we can't determine this during decoration - see, e.g.,
                http://stackoverflow.com/questions/2366713/ )
            outputs (list or None): NOT IMPLEMENTED. optional; a list of output files to retrieve from remote worker.
            cfg (dict): NOT IMPLEMENTED. Configuration to use on remote worker.
        """
        # If the wrapped function doesn't ALWAYS run remotely,
        # then it has to return its standard return value
        self.outputs = outputs

        cfg = if_not_none(cfg, config.cfg)

        self.remote = remote
        self.display = display
        self.sendsource = sendsource
        self.image = if_not_none(image, cfg.default_image)
        self.engine = if_not_none(engine, cfg.default_engine)
        self.job_config = if_not_none(config.cfg, {})
        self.jobname = jobname
        self.is_imethod = is_imethod

    def __call__(self, func):
        """
        This gets called with the function we wish to wrap
        """

        assert callable(func)

        if self.jobname is None:
            self.jobname = func.__name__

        @utils.args_from(func)  # TODO: should inject 'wait', deal with functools.wraps
        def wrapper(*args, **kwargs):
            """This documentation should be replaced with `func`'s by the decorator"""

            f = func  # because we reassign instance methods
            if not wrapper.remote:
                return f(*args, **kwargs)

            # TODO: this doesn't work because the decorator gives wrapper f's signature
            wait = kwargs.get('wait', True)

            # Bind methods to their objects
            if self.is_imethod:
                # We can't call this function like normal, because the decorators can't identify instance methods
                # Instead, we'll create another bound copy of the instancemethod (probably only need to do this once)
                fn_self = args[0]
                f = types.MethodType(f, fn_self, fn_self.__class__)
                args = args[1:]

            # Submit job to remote engine
            python_call = bpy.PythonCall(f, args, kwargs)
            job = bpy.PythonJob(self.engine,
                                self.image,
                                python_call,
                                name=self.jobname,
                                sendsource=self.sendsource)
            if self.display:
                mdt.logs.display(job, title=f.__name__)

            if wait:
                job.wait()
                return job.finish()
            else:
                return job

        wrapper.remote = self.remote
        return wrapper


runsremotely = RunsRemotely  # because decorators should be lower case

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