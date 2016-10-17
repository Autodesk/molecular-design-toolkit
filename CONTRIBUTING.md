# Contributing to Molecular Design Toolkit


<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Contributing to Molecular Design Toolkit](#contributing-to-molecular-design-toolkit)
    - [Tips and guidelines](#tips-and-guidelines)
        - [What can I contribute?](#what-can-i-contribute)
        - [Pull requests are always welcome](#pull-requests-are-always-welcome)
        - [Design and cleanup proposals](#design-and-cleanup-proposals)
    - [Submission Guidelines](#submission-guidelines)
        - [Project Roles](#project-roles)
        - [Timing](#timing)
        - [Issues](#issues)
        - [Pull Requests](#pull-requests)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

<!-- to generate: npm install doctoc: doctoc --gitlab --maxlevel 3 CONTRIBUTING.md-->



## Tips and guidelines

### What can I contribute?
Contributions to this project are encouraged! Email the maintainers at `moldesign_maintainers@autodesk.com` to become a contributor.

If you're interested in getting started, here are some general contributions that are always welcome. We also maintain a [wishlist of specific ideas on the wiki](https://github.com/Autodesk/molecular-design-toolkit/wiki/Contribution-ideas).

**Tests**: This is one of the easiest ways to get started - see also `moldesign/tests/README.md`

 * **Contribute unit tests** - test that MDT's functions work correctly in a variety of situations, and that they report errors when appropriate 
 * **Validate methods** - methods need to be tested against known results for accuracy; for instance, we need to check that a Hartree-Fock calculation gives the same results as Gaussian, GAMESS, and MolDesign

**Examples**: See also `moldesign/notebooks`

 * **Best practices** - put together template notebooks that will help users get started with a particular workflow, from modeling proteins from crystal structures to exploring electronic structure
 * **Interesting use cases** - Contribute a notebook that shows off a cool bit of science or design.

**Bug fixes:** Found a typo in the code? Found that a function fails under certain conditions? Know how to fix it? Great! Go for it. Please do [open an issue](https://github.com/autodesk/molecular-design-toolkit/issues) so that we know you're working on it, and submit a pull request when you're ready.

**Features:** The chemical modeling universe is vast, and we want toolkit users to have access to a lot of it. Whether you want free energy perturbation, or Boyes' localization, or 3D structural alignment - we want it too! As always, please [open an issue](https://github.com/autodesk/molecular-design-toolkit/issues) so that we know what you're working on.


**Whatever:** There's ALWAYS something to do, whether supporting other languages (e.g., Spanish or Bahasa Indonesia, not Fortran or C++); improving 3D viewer performance; improving documentation; adding Python 3.X support; or integrating with other IDE technologies.

### Pull requests are always welcome

All PRs should be documented as [GitHub issues](https://github.com/autodesk/molecular-design-toolkit/issues), ideally BEFORE you start working on them.

### Design and cleanup proposals

Good API design is at the heart of this project, and you don't need to do any programming to help with this! For example:

 * You could describe how a user will run a Hamiltonian replica exchange calculation (should it be a class or a function? What are the method names? How does the user specify the temperatures?).
 * You can also propose redesigns for existing features - maybe you think `mdt.add_hydrogens(mol)` should be renamed to `mol.add_hydrogens()`, or you want to propose a better way to access trajectory data.

To get started, as always: [open an issue](https://github.com/autodesk/molecular-design-toolkit/issues). For information on making these types of
contributions, see [the development guide](DEVELOPMENT.md).



## Submission Guidelines

### Maintainers
Maintainers are responsible for responding to pull requests and issues, as well as guiding the direction of the project.

Aaron Virshup - Lead developer and maintainer<br>
Dion Amago - Maintainer<br>
Malte Tinnus - Maintainer

If you've established yourself as an impactful contributor for the project, and are willing take on the extra work, we'd love to have your help maintaining it! Email the maintainers list at `moldesign_maintainers@autodesk.com` for details.

### Timing

We will attempt to address all issues and pull requests within one week. It may a bit longer before pull requests are actually merged, as they must be inspected and tested. 

### Issues

If MDT isn't working like you expect, please open a new issue! We appreciate any effort you can make to avoid reporting duplicate issues, but please err on the side of reporting the bug if you're not sure.

Providing the following information will increase the chances of your issue being dealt with quickly:

* **Overview of the Issue** - Please describe the issue, and include any relevant exception messages or screenshots.
* **Environment** - Include the relevant output of `pip freeze` as well as your system and python version info.
* **Help us reproduce the issue** - Please include code that will help us reproduce the issue. For complex situations, attach a notebook file.
* **Related Issues** - Please link to other issues in this project (or even other projects) that appear to be related 

### Pull Requests

Before you submit your pull request consider the following guidelines:


* Search GitHub for an open or closed Pull Request that relates to your submission. You don't want to duplicate effort.
* Make your changes in a new git branch:

     ```shell
     git checkout -b my-fix-branch [working-branch-name]
     ```

* Create your patch.
* Commit your changes using a descriptive commit message.

     ```shell
     git commit -a
     ```
  Note: the optional commit `-a` command line option will automatically "add" and "rm" edited files.

* Push your branch to GitHub:

    ```shell
    git push origin my-fix-branch
    ```

* In GitHub, send a pull request to `molecular-design-toolkit:dev`
* Before any request is merged, you'll need to agree to the contribution boilerplate. Email us at `moldesign_maintainers@autodesk.com` for details. 
