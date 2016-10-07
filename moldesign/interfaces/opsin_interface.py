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

def name_to_smiles(name,
                   image='opsin',
                   engine=None):

    command = 'opsin -osmi input.txt output.txt'

    # TODO: this boilerplate has to go
    engine = utils.if_not_none(engine, mdt.compute.get_engine())
    imagename = mdt.compute.get_image_path(image)
    job = engine.launch(imagename,
                          command,
                          inputs={'input.txt': name + '\n'},
                          name="opsin, %s" % name)
    mdt.uibase.display_log(job.get_display_object(), "opsin, %s"%name)
    job.wait()
    smistring = job.get_output('output.txt').read().strip()
    if not smistring:
        raise ValueError('Could not parse chemical name "%s"' % name)
    else:
        return smistring
