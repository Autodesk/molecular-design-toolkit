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

import types

from pyccc import python as bpy

import moldesign as mdt
from moldesign import utils, uibase
from . import configuration

class RunsRemotely(object):
    def __init__(self, enable=True,
                 display=True,
                 jobname=None,
                 sendsource=False,
                 engine=None,
                 image=None,
                 is_imethod=False):
        """Function decorator to run a python function remotely.

        Note:
         This ONLY works for pure functions - where you're interested in the
           return value only. Side effects won't be visible to the user.

        Args:
            enable (bool): If True, run this job using a compute engine
            display (bool): Create a jupyter logging display for the remote job
                (default: True in Jupyter notebooks, False otherwise)
            jobname (str): Name metadata - defaults to the __name__ of the function
            sendsource (bool): if False (default), call this function directly on the remote worker;
               if True, send the function's source code (for debugging, mostly)
            engine (pyccc.engine.EngineBase): engine to send the job to (default:
                moldesign.compute.get_engine())
            image (str): name of the docker image (including registry, repository, and tags)
                (default: moldesign.config.default_python_image)
            is_imethod (bool): This is an instancemethod
               Note: we can't determine this at import-time without going to great lengths ...
                    - see, e.g., http://stackoverflow.com/questions/2366713/ )
        """
        self.enabled = enable
        self.display = display
        self.sendsource = sendsource
        self.image = image
        self.engine = engine
        self.jobname = jobname
        self.is_imethod = is_imethod

    def __call__(self, func):
        """
        This gets called with the function we wish to wrap
        """
        assert callable(func)

        if self.jobname is None:
            self.jobname = func.__name__

        if func.__name__ == 'wrapper': assert False

        @utils.args_from(func,
                         wraps=True,
                         inject_kwargs={'wait': True})
        def wrapper(*args, **kwargs):
            """ Wraps a python function so that it will be executed remotely using a compute engine

            Note:
                At runtime, this documentation should be replaced with that of the wrapped function
            """
            # If the wrapper is not enabled, just run the wrapped function as normal.
            f = func  # keeps a reference to the original function in this closure
            if not wrapper.enabled:
                return f(*args, **kwargs)

            wait = kwargs.get('wait', True)

            # Bind instance methods to their objects
            if self.is_imethod:
                f, args = _bind_instance_method(f, args)

            # Submit job to remote engine
            python_call = bpy.PythonCall(f, *args, **kwargs)
            image = utils.if_not_none(self.image, configuration.config.default_python_image)
            engine = utils.if_not_none(self.engine, mdt.compute.get_engine())
            job = bpy.PythonJob(engine,
                                image,
                                python_call,
                                name=self.jobname,
                                sendsource=self.sendsource)

            if self.display:
                uibase.display_log(job.get_display_object(), title=f.__name__)

            if wait:
                job.wait()
                return job.result
            else:
                return job

        wrapper.__name__ = func.__name__
        wrapper.enabled = self.enabled
        return wrapper

runsremotely = RunsRemotely  # because decorators should be lower case


def _bind_instance_method(f, args):
    # We can't call this function like normal, because the decorators can't identify
    # instance methods. Instead, we'll create another bound copy of the instancemethod (probably
    # only need to do this once)
    fn_self = args[0]
    f = types.MethodType(f, fn_self, fn_self.__class__)
    args = args[1:]
    return f, args

