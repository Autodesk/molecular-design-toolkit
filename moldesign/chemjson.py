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

import json
import numpy


# TODO: register deserializers (with decorators?)
class JsonEncoder(json.JSONEncoder):
    def default(self, obj):
        if hasattr(obj, 'to_json'):
            obj = obj.to_json()

        return obj


def writer(obj, fileobj):
    """ Serialize an object to JSON.

    This is quite general - the object must be valid JSON or must have a "to_json" method.

    Note:
        Usually not called directly - usually called from mdt.write()

    Args:
        obj (object): object to serialize; must be valid JSON or have a to_json method
        fileobj (filelike): File-like object to serialize to
    """
    json.dump(obj, fileobj, cls=JsonEncoder)


def reader(fileobj):
    pass


def jsonify(obj, attrnames):
    js = {}
    for item in attrnames:
        attr = getattr(obj, item, None)
        if attr is None: continue
        js[item] = attr

    return js