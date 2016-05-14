# The Molecular Design Toolkit

Molecular modeling without the pain - a Python 2.7 library offering integrated simulation, visualization, analysis, and cloud computing. 

The toolkit aims to lower the barriers between you and your science by integrating mature, open source simulation packages with a readable abstract API, Jupyter notebook visualization, and native cloud computing.

**Install it**: `pip install moldesign`
**Launch an example notebook**: `python -m moldesign intro`

## Code Example

You'll almost always import the package and its units module:
<pre><code>import moldesign as mdt
from moldesign import units as u
</code></pre>

Download a protein from the PDB and visualize it in 3D (in a notebook):
<pre><code>
protease = mdt.from_pdb('3AID')
protease.draw()
</code></pre>

Create a small molecule and relax its geometry:
<pre><code>mol = mdt.from_name('bipyridine')
mol.set_energy_model(mdt.models.RHF(basis='STO-3G'))
min_trajectory = mol.minimize(nsteps=20)
min_trajectory.draw_orbitals()
</code></pre>

For in-depth examples, see the built-in example notebooks (run `python -m moldesign intro` to launch).


## Documentation

API documentation link
User documentation link
Forums


## Tests

`moldesign/tests` 

## Contributors

CONTRIBUTING.md

## License

Copyright 2015 Autodesk Inc.

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.