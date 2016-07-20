# Creating new MDT project releases

These instructions generally apply to [Autodesk/molecular-design-toolkit](https://github.com/Autodesk/molecular-design-toolkit),
[Autodesk/py-cloud-compute-cannon](https://github.com/Autodesk/py-cloud-compute-cannon), and [Autodesk/notebook-molecular-visualization](https://github.com/Autodesk/notebook-molecular-visualization).

This is only for project maintainers - everyone else should use PRs.

Everything here is HIGHLY PRELIMINARY!!! This will change dramatically once Jenkins or Travis is up.

1. Push changes or merge PR into dev branch and check out a clean copy.
1. Test it.
1. Run `check-manifest` at project root. (Repeat steps 1-3 until this checks out OK)
1. [MDT only] - update/push dockerfiles and image URLs in `moldesign.compute.configuration`
2. Pull-request into master
3. Tag master branch with new release number
4. In repo root. `python setup.py register -r pypi`
5. In repo root, `python setup.py sdist upload -r pypi`
6. [MDT only] Run `cd docs; make html` and push contents of `_build/html` to the gh-pages branch.