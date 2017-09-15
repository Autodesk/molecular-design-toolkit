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
import types
import future.utils
from pyccc import python as bpy

import moldesign as mdt
from moldesign import utils
from . import configuration, run_job
from ..helpers import display_log


class RpcWrapper(object):
    """ A wrapper that lets to transparently execute python functions in remote
    environments - usually in docker containers.

    These wrappers are built to allow a lot of run-time flexibility based on the description
    of the package (``self.pkg``) that's being called.

    Note:
     This ONLY works for pure functions - where you're interested in the
       return value only. Side effects - including any object state - will be discarded.

    Args:
        pkg (mdt.compute.packages.InterfacedPackage): package to run this command with
        display (bool): Create a jupyter logging display for the remote job
            (default: True in Jupyter notebooks, False otherwise)
        jobname (str): Name metadata - defaults to the __name__ of the function
        sendsource (bool): if False (default), call this function directly on the remote worker;
           if True, send the function's source code (for debugging, mostly)
        persist_refs (bool): Persist python object references across the RPC roundtrip
        is_imethod (bool): This is an instancemethod
           Note: we can't determine this at import-time without going to great lengths ...
                - see, e.g., http://stackoverflow.com/questions/2366713/ )
    """
    def __init__(self, pkg,
                 display=True,
                 jobname=None,
                 sendsource=False,
                 is_imethod=False,
                 persist_refs=False):
        self.pkg = pkg
        self.display = display
        self.sendsource = sendsource
        self.jobname = jobname
        self.is_imethod = is_imethod
        self.persist_refs = persist_refs

    def __call__(self, func):
        """
        This gets called with the function we wish to wrap
        """
        from .compute import get_image_path
        assert callable(func)

        if self.jobname is None:
            self.jobname = func.__name__

        assert func.__name__ != 'wrapper'  # who wraps the wrappers?

        @utils.args_from(func,
                         wraps=True,
                         inject_kwargs={'wait': True})
        def wrapper(*args, **kwargs):
            """ Wraps a python function so that it will be executed remotely using a compute engine

            Note:
                At runtime, this documentation should be replaced with that of the wrapped function
            """
            f = func  # keeps a reference to the original function in this closure
            wait = kwargs.get('wait', True)

            if wait and not self.pkg.force_remote:
                return f(*args, **kwargs)

            # Bind instance methods to their objects
            if self.is_imethod:
                f, args = _bind_instance_method(f, args)

            # Submit job to remote engine
            python_call = bpy.PythonCall(f, *args, **kwargs)

            engine = utils.if_not_none(self.pkg.engine, mdt.compute.get_engine())
            job = bpy.PythonJob(engine=engine,
                                image=self.pkg.get_docker_image_path(),
                                command=python_call,
                                name=self.jobname,
                                sendsource=self.sendsource,
                                interpreter='python',  # always run in image's native interpreter
                                persist_references=self.persist_refs,
                                submit=False)

            return run_job(job, wait=wait, _return_result=True)

        wrapper.__name__ = func.__name__
        wrapper.__wrapped__ = func
        return wrapper


def _bind_instance_method(f, args):
    # We can't call this function like normal, because the decorators can't identify
    # instance methods. Instead, we'll create another bound copy of the instancemethod (probably
    # only need to do this once)
    fn_self = args[0]
    if future.utils.PY2 == 2:
        f = types.MethodType(f, fn_self, fn_self.__class__)
    else:
        f = types.MethodType(f, fn_self)
    args = args[1:]
    return f, args

