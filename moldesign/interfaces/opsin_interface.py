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
import pyccc

import moldesign as mdt
from moldesign import utils

IMAGE = 'opsin'


@utils.kwargs_from(mdt.compute.run_job)
def name_to_smiles(name,
                   **kwargs):

    command = 'opsin -osmi input.txt output.txt'

    def finish_job(job):
        smistring = job.get_output('output.txt').read().strip()
        if not smistring:
            raise ValueError('Could not parse chemical name "%s"' % name)
        else:
            return smistring

    job = pyccc.Job(image=mdt.compute.get_image_path(IMAGE),
                    command=command,
                    name="opsin, %s" % name,
                    inputs={'input.txt': name + '\n'},
                    when_finished=finish_job)

    return mdt.compute.run_job(job, _return_result=True, **kwargs)
