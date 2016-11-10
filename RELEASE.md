# Whenever commiting changes anywhere

Make SURE that you've run `set_filters.sh` at the project base. This will make sure that you don't accidentally commit any `ipynb` output fields into the repository. You can check to make sure the filters are working by running `git diff` on any notebook file that has output in it: all `output` and `metadata` fields should remain blank.


# Maintainers: Accepting a PR

Work can be merged into `master` or a feature branch, as appropriate. Don't merge broken code
into master, but it doesn't need to be totally production-ready either: only releases (below)
are considered "stable".

1. Review the code.
1. Make sure that there's appropriate functional tests.
1. Check that the travis build is at least running all the way to the end. The tests don't *necessarily* need to pass, but you need to undertand why  what's passing and what's not.


# Maintainers: Creating a new release

1. Decide on the new version number (using [semantic versioning](http://semver.org/)). For our purposes here, we'll pretend it's `0.9.3`.
1. Tag the relevant commit (the build must be passing) with a release candidate version number, e.g., `0.9.3rc1`. Note that this commit must be passing the full nightly test battery in all test environments (still WIP 11/2/16)
1. After travis finishes building all deployment artifacts (still WIP 11/2/16), manually test all examples and tutorials from a docker container by:
  A. Running `docker run -it -p 8888:8888 moldesign_notebook`, then
  B. Testing the notebooks at http://localhost:8888
1. If something isn't working right, keep working, and keep tagging release candidates with `0.9.3.rc2`, `0.9.3rc3`, ..., until it works
1. If it IS working, tag THE SAME COMMIT YOU JUST TESTED (the one already tagged as a release candidate) with its final version string - that's `0.9.3` here.


# Maintainers: updating docs

Documentation is NOT coupled to the package releases; docs tend to get updated continuously.

1. In the `master` branch, update the version numbers in `docs/conf.py`
1. Run `cd docs; make clean; make html`. 
1. In a separate directory, check out a fresh copy of the repo and run `git checkout gh-pages`
1. Copy the contents of `[master branch]/docs/_build/html` into the root of the `gh-pages` branch.
1. Commit your changes to the `gh-pages` branch and push them back to GitHub.
