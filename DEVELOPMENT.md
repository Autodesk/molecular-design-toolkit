# DEVELOPING MDT

### Setting up a dev environment
(still under construction)

### Install prequisites (first time only)
You need to install docker, and an environment manager for Python 3 (Miniconda 3). Here's one way to do that:
1. Install docker: [link]
2. Install pyenv and pyenv-venv: [link]
3. Install miniconda3 by running: `pyenv install miniconda3-latest`
4. Switch to miniconda environment by running: `pyenv shell miniconda3-latest`

### Set up your environment (first time only)
1. Get MDT: `git clone http://github.com/Autodesk/molecular-design-toolkit`
1. `cd molecular-design-toolkit`
1. Create conda environment (optional but recommended) by running: [command to create conda env]
2. Activate the environment: `pyenv activate [environment name???]`
1. Install dev dependencies: `pip install -r requirements.txt DockerMakefiles/requirements.txt deployment/requirements.txt`
2. Set up for local dev mode (this tells MDT to use your local docker containers):
```bash
mkdir ~/.moldesign
echo "devmode: true" > ~/.moldesign/moldesign.yml
```
8. Install MDT in "development mode": 
```
pip install -e molecular-design-toolkit
```

### To activate environment (in any new shell)
1. Run `pyenv activate [environment name???]`

### To rebuild docker images (first time and after changes that affect dockerized code)
5. Build development versions of all docker images:
```bash
cd DockerMakefiles
docker-make --all --tag dev
```

### To run tests
```bash
cd molecular-design-toolkit/moldesign/_tests
py.test -n [number of concurrent tests]
```

See [the testing README](moldesign/_tests/README.md) for more details.


 
### Code style
1. Functions and variables should be `lowercase_with_underscores`. Class names and constructors should be `CapitalizedCamelCase`.
1. The user-exposed API should be clean, [PEP 8](https://www.python.org/dev/peps/pep-0008/) code.
1. Internally, readability and functionality are more important than consistency - that said,  [Google's style guide](https://google.github.io/styleguide/pyguide.html), along with [PEP 8](https://www.python.org/dev/peps/pep-0008/), is strongly encouraged.


# Contributing

### Who should contribute?
Anyone with a molecular modeling workflow that they want to enable or share. Experience and research-level knowledge of the field is an important asset! In contrast, limited programming experience *is definitely not* barrier to contributing - we can help you with that! Please ask for help getting started in our forums [link].

### Scope: What goes into MDT?
Established techniques and general simulation tools that will be useful for **3-dimensional biomolecular modeling**. MDT aims to enable scientists to easily build new simulation techniques and workflows, but new, immature techniques, or those with limited applicability outside of a particular system should be implemented as separate projects that *use* MDT, not *part of* MDT.
###### Could (and should!) be implemented in MDT:
 * Physical simulation and modelilng: Lambda dynamics; homology modelling; surface hopping; RPMD; metadynamics; markov state models; a library of common 3D structures (such as amino acids, carbon nanotubes, small molecules, etc.)
 * Visualization and UI: transitions between different views; interactive structure building and editing; ray-traced rendering; movie exports
 
###### Should implemented as a separate project:
 * Computational techniques: fluid dynamics solver (not useful at the atomic level), biological network models (no clear connection to 3D structures); machine-learning based quantum chemistry (immature, untested)
 * Visualization and UI: visualizations for specific systems (not generally applicable); 


# Development guidelines

### Whenever commiting changes anywhere

Make SURE that you've run `nb-output-filter.sh` at the project base. You only need to do this once (per copy of the repository). This will make sure that you don't accidentally commit any `ipynb` output fields into the repository. You can check to make sure the filters are working by running `git diff` on any notebook file that has output in it: all `output` and `metadata` fields should remain blank.


### Maintainers: Accepting a PR

Work can be merged into `master` or a feature branch, as appropriate. Don't merge broken code
into master, but it doesn't need to be totally production-ready either: only releases (below)
are considered "stable".

1. Review the code.
1. Make sure that there's appropriate functional tests.
1. Check that the travis build is at least running all the way to the end. The tests don't *necessarily* need to pass, but you need to undertand why  what's passing and what's not.


### Releases

1. Decide on the new version number (see below). For our purposes here, we'll pretend it's `0.9.3`.
1. Tag the relevant commit (the build must be passing) with a release candidate version number, e.g., `0.9.3rc1`.
1. Codeship will automatically deploy the updated release to PyPI and DockerHub
1. Manually test the example notebooks against this pre-release version.
1. If succesful, tag the relevant commit with the official release version `0.9.3`

### Versioning
For now, we're using a subset [PEP 440](https://www.python.org/dev/peps/pep-0440/):
1. Every release should be of the form MAJOR.MINOR.PATCH, e.g. `0.1.2`
2. Pre-releases should be numbered consecutively, and may be alpha, beta, or "release candidate", e.g. `1.0.1rc3` or `0.5.3a1`
3. Our deployment infrastructure uses this regular expression to accept version strings:
`^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)((a|rc|b)(0|[1-9]\d*))?$`

### Maintainers: updating the documentation

Documentation is NOT coupled to the package releases; docs tend to get updated continuously.

1. In the `master` branch, update the version numbers in `docs/conf.py`
1. Run `cd docs; make clean; make html`. 
1. In a separate directory, check out a fresh copy of the repo and run `git checkout gh-pages`
1. Copy the contents of `[master branch]/docs/_build/html` into the root of the `gh-pages` branch.
1. Commit your changes to the `gh-pages` branch and push them back to GitHub.
