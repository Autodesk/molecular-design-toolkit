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
import os
from functools import wraps

import collections
import funcsigs

from .utils import if_not_none
from .docparsers import GoogleDocArgumentInjector

def args_from(original_function,
              only=None,
              allexcept=None,
              inject_kwargs=None,
              wraps=None,
              append_docstring_description=False,
              transfer_docstring_args=False):
    """
    Decorator to transfer call signatures - helps to hide ugly *args and **kwargs in delegated calls

    Note:
        Docstring modification is not implemented. It will be intended to work with google-style
        docstrings only.

    Args:
        original_function (callable): the function to take the call signature from
        only (List[str]): only transfer these arguments (incompatible with `allexcept`)
        wraps (bool): Transfer documentation and attributes from original_function to
            decorated_function, using functools.wraps (default: True if call signature is
            unchanged, False otherwise)
        allexcept (List[str]): transfer all except these arguments (incompatible with `only`)
        inject_kwargs (dict): Inject new kwargs into the call signature
            (of the form `{argname: defaultvalue}`)
        transfer_docstring_args (bool): Update the 'Args:' section of the docstring with the
            previous function's docs
        append_docstring_description (bool): Update all sections of the docstring (other than args)

    Returns:
        Decorator function
    """
    # TODO - deal with docstrings - how does wraps interact with append/transfer_docstring?
    # NEWFEATURE - verify arguments?

    if append_docstring_description or transfer_docstring_args:
        raise NotImplementedError()

    if only and allexcept:
        raise ValueError('Error in keyword arguments - '
                         'pass *either* "only" or "allexcept", not both')

    if hasattr(original_function, '__signature__'):
        sig = original_function.__signature__.replace()
    else:
        sig = funcsigs.signature(original_function)

    # Modify the call signature if necessary
    if only or allexcept or inject_kwargs:
        wraps = if_not_none(wraps, False)
        newparams = []
        if only:
            for param in only:
                newparams.append(sig.parameters[param])
        elif allexcept:
            for name, param in sig.parameters.iteritems():
                if name not in allexcept:
                    newparams.append(param)
        else:
            newparams = sig.parameters.values()
        if inject_kwargs:
            for name, default in inject_kwargs.iteritems():
                newp = funcsigs.Parameter(name, funcsigs.Parameter.POSITIONAL_OR_KEYWORD,
                                          default=default)
                newparams.append(newp)

        sig = sig.replace(parameters=newparams)

    else:
        wraps = if_not_none(wraps, True)

    def decorator(f):
        """Modify f's call signature (using the `__signature__` attribute)"""
        fname = f.__name__
        if wraps:
            f = functools.wraps(original_function)(f)
            f.__name__ = fname  # revert name change
        f.__signature__ = sig

        # Modify docstring for documentation only
        if os.environ.get('SPHINX_IS_BUILDING_DOCS', ""):
            sigstring = '%s%s\n' % (f.__name__, sig)
            if hasattr(f, '__doc__') and f.__doc__ is not None:
                f.__doc__ = sigstring + f.__doc__
            else:
                f.__doc__ = sigstring
        return f

    return decorator


def kwargs_from(reference_function, mod_docs=True):
    """ Replaces ``**kwargs`` in a call signature with keyword arguments from another function.

    Args:
        reference_function (function): function to get kwargs from
        mod_docs (bool): whether to modify the decorated function's docstring

    Note:
        ``mod_docs`` works ONLY for google-style docstrings
    """
    refsig = funcsigs.signature(reference_function)

    kwparams = []
    for name, param in refsig.parameters.iteritems():
        if param.default != param.empty or param.kind in (param.VAR_KEYWORD, param.KEYWORD_ONLY):
            kwparams.append(param)

    if mod_docs:
        refdocs = GoogleDocArgumentInjector(reference_function.__doc__)

    def decorator(f):
        sig = funcsigs.signature(f)

        fparams = []
        found_varkeyword = None

        for name, param in sig.parameters.iteritems():
            if param.kind == param.VAR_KEYWORD:
                fparams.extend(kwparams)
                found_varkeyword = name
            else:
                fparams.append(param)

        if not found_varkeyword:
            raise TypeError("Function has no **kwargs wildcard.")

        f.__signature__ = sig.replace(parameters=fparams)

        if mod_docs:
            docs = GoogleDocArgumentInjector(f.__doc__)
            new_args = collections.OrderedDict()
            for argname, doc in docs.args.iteritems():
                if argname == found_varkeyword:
                    for name in kwparams:
                        new_args[name] = refdocs.args[name]
                else:
                    new_args[argname] = doc

            f.__doc__ = docs.new_docstring()
        return f

    return decorator


class DocInherit(object):
    """
    Allows methods to inherit docstrings from their superclasses
    FROM http://code.activestate.com/recipes/576862/
    """
    def __init__(self, mthd):
        self.mthd = mthd
        self.name = mthd.__name__

    def __get__(self, obj, cls):
        if obj:
            return self.get_with_inst(obj, cls)
        else:
            return self.get_no_inst(cls)

    def get_with_inst(self, obj, cls):

        overridden = getattr(super(cls, obj), self.name, None)

        @wraps(self.mthd, assigned=('__name__','__module__'))
        def f(*args, **kwargs):
            return self.mthd(obj, *args, **kwargs)

        return self.use_parent_doc(f, overridden)

    def get_no_inst(self, cls):

        for parent in cls.__mro__[1:]:
            overridden = getattr(parent, self.name, None)
            if overridden: break

        @wraps(self.mthd, assigned=('__name__','__module__'))
        def f(*args, **kwargs):
            return self.mthd(*args, **kwargs)

        return self.use_parent_doc(f, overridden)

    def use_parent_doc(self, func, source):
        if source is None:
            raise NameError, ("Can't find '%s' in parents"%self.name)
        func.__doc__ = source.__doc__
        return func

#idiomatic decorator name
doc_inherit = DocInherit