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
class Categorizer(dict):
    """
    Create a dict of lists from an iterable, with dict keys given by keyfn
    """

    def __init__(self, keyfn, iterable):
        super(Categorizer, self).__init__()
        self.keyfn = keyfn
        for item in iterable:
            self.add(item)

    def add(self, item):
        key = self.keyfn(item)
        if key not in self:
            self[key] = []
        self[key].append(item)


class DotDict(dict):
    """Dict with items accessible as attributes"""
    def __getstate__(self):
        retval = dict(__dict__=self.__dict__.copy(),
                      items=self.items())
        return retval

    def __setstate__(self, d):
        self.__dict__.update(d['__dict__'])
        self.update(d['items'])

    def __getattr__(self, item):
        try:
            return self[item]
        except KeyError:
            raise AttributeError("'%s' object has no attribute '%s'"
                                 % (self.__class__.__name__, item))

    def __setattr__(self, item, val):
        self[item] = val

    def __dir__(self):
        return dir(self.__class__) + self.keys()


class Alias(object):
    """
    Descriptor that calls a child object's method.
    e.g.
    >>> class A(object):
    >>>     childkeys = Alias('child.keys')
    >>>     child = dict()
    >>>
    >>> a = A()
    >>> a.child['key'] = 'value'
    >>> a.childkeys() #calls a.child.keys(), returns ['key']
    ['key']
    """
    def __init__(self, objmethod):
        objname, methodname = objmethod.split('.')
        self.objname = objname
        self.methodname = methodname

    def __get__(self, instance, owner):
        proxied = getattr(instance, self.objname)
        return getattr(proxied,self.methodname)


class Synonym(object):
    """ An attribute (class or intance) that is just a synonym for another.
    """
    def __init__(self, name):
        self.name = name

    def __get__(self, instance, owner):
        return getattr(instance, self.name)

    def __set__(self, instance, value):
        setattr(instance, self.name, value)


class DictLike(object):
    """
    This just wraps normal dicts so that other classes don't have to inherit from a built-in class,
    which apparently breaks pickle quite frequently.
    """
    def __getstate__(self):
        retval = dict(__dict__=self.__dict__.copy())
        return retval

    def __setstate__(self, d):
        self.__dict__.update(d['__dict__'])

    def __init__(self, **kwargs):
        self.children = {}
        self.children.update(kwargs)

    def __getattr__(self, k):
        return getattr(self.children, k)

    __setitem__ = Alias('children.__setitem__')
    __getitem__ = Alias('children.__getitem__')

    def __len__(self):
        return len(self.children)


class Attribute(object):
    """For overriding a property in a superclass - turns the attribute back
    into a normal instance attribute"""
    def __init__(self, name):
        self.name = name

    def __get__(self, instance, cls):
        return getattr(instance, self.name)

    def __set__(self, instance, value):
        return setattr(instance, self.name, value)