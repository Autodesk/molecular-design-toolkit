# Whenever commiting changes anywhere

Make SURE that you've run `set_filters.sh` at the project base. This will make sure that you don't accidentally commit any `ipynb` output fields into the repository. You can check to make sure the filters are working by running `git diff` on any notebook file that has output in it: all `output` and `metadata` fields should remain blank.

# Creating a Pull Request

These instructions generally apply to [Autodesk/molecular-design-toolkit](https://github.com/Autodesk/molecular-design-toolkit),
[Autodesk/py-cloud-compute-cannon](https://github.com/Autodesk/py-cloud-compute-cannon), and [Autodesk/notebook-molecular-visualization](https://github.com/Autodesk/notebook-molecular-visualization).

1. Group all changes in a single branch.
1. Run tests in `moldesign/_tests` with `py.test -n NCORES` where `NCORES` is the number of jobs to run simultanously.
2. Create a pull-request for the appropriate branch in Autodesk/molecular-design-toolkit


# Maintainers: Accepting a PR

Work can be merged into `master` or a feature branch, as appropriate. Don't merge broken code
into master, but it doesn't need to be totally production-ready either: only releases (below)
are considered "stable".

1. Review the code.
1. Make sure that there's appropriate functional tests.
1. Run all tests. They don't all *necessarily* need to pass, but you need to undertand why  what's passing and what's not.
1. Run any example notebooks that might be affected by the PR.


# Maintainers: Creating a new release

Decide on the new version number (using [semantic versioning](http://semver.org/)).

Everything here is HIGHLY PRELIMINARY!!! This will be a lot easier once Travis/Jenkins is up.

1. Make sure all changes have been PR'd into master
1. Check out a clean copy of master
1. Run `check-manifest` at project root to make sure the distribution will include the necessary files.
1. Increment the `default_version_tag` field in `moldesign.compute.config`, and commit the change
1. Tag your local branch with the release number: `git tag [version]` (note: delete this tag with `git tag rm [version]` if you need to abort the release)
4. Build new docker images: `cd docker_images; ./dockermake.py --all --repo docker.io/Autodesk/moldesign: --tag [version]`
1. Confirm that all tests are passing *with and without* locally installed dependencies (pyscf, openbabel, etc)
1. Confirm that all tutorial and example notebooks run without errors.

If this is all succesfull - you're ready to make it public.

1. Push docker images to cloud: `cd docker_images; ./dockermake.py --all --push --repo docker.io/Autodesk/moldesign: --tag [version]`
4. `python setup.py register -r pypi`
5. `python setup.py sdist upload -r pypi`
1. `git push origin master --tags`

The final step is the point of no return - you'll need to prep a new release if you discover a problem. Before that, you can undo
what you've done by deleting your release off of PyPI and DockerHub.



# Maintainers: updating docs

Documentation is NOT coupled to the package releases; docs tend to get updated continuously.

1. In the `master` branch, update the version numbers in `docs/conf.py`
1. Run `cd docs; make clean; make html`. 
1. In a separate directory, check out a fresh copy of the repo and run `git checkout gh-pages`
1. Copy the contents of `[master branch]/docs/_build/html` into the root of the `gh-pages` branch.
1. Commit your changes to the `gh-pages` branch and push them back to GitHub.
