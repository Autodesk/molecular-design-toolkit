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
from ..utils import exports
import functools
import imp

try:
    imp.find_module('nbmolviz')
except ImportError:
    nbmolviz_installed = False
else:
    nbmolviz_installed = True

__init = False


def _get_nbmethod(name):
    # don't import nbmolviz methods until a method is actually called
    from nbmolviz import methods as nbmethods
    module = nbmethods
    for item in name.split('.'):
        module = getattr(module, item)
    return module


@exports
class WidgetMethod(object):
    def __init__(self, name):
        self.name = name
        self.method = None

    def __get__(self, instance, owner):
        if self.method is None:
            self.method = _get_nbmethod(self.name)
        return functools.partial(self.method, instance)


def not_installed_method(*args, **kwargs):
    raise ImportError(
            "To use MDT's graphical user interfaces in a Noteobook, please install "
            "the nbmolviz library: `pip install nbmolviz`")


