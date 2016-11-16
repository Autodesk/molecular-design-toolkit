**************
Installing MDT
**************

.. contents::
  :depth: 2



Basic installation
==================
MDT requires Python 2.7. Type

.. code-block:: bash

    python --version

to make sure that you have Python 2.7.X.  If not, you'll want to use your `system's package manager. <#system-specific-installation>`_ If have Python 2.7 but not ``pip`` (the python package manager), you can install it by running

.. code-block:: bash

    easy_install pip

At this point, you should have everything you need to install MDT.

.. code-block:: bash

    pip install moldesign

**Note:** Depending on their specific python installation, some users may need to run the installation as root, i.e. ``sudo pip install moldesign``



Updating
^^^^^^^^

To update to the most recent version of the toolkit, type

.. code-block:: bash

    pip install -U moldesign

However, note that this will update ``moldesign`` and *all of its dependencies* (including `numpy`, `jupyter`, `setuptools`, etc.)  to their latest versions, which you may or may not want. To update only the ``moldesign`` suite, please run

.. code-block:: bash

    pip install --no-deps -U moldesign pyccc nbmolviz


Common problems
^^^^^^^^^^^^^^^

**Permissions**
Depending on how python is installed on your system, you may need to run this installation as root, e.g. ``sudo pip install moldesign``.

**MacOS default Python**
We've encountered issues trying to install Jupyter with MacOS's built-in python distribution. We highly recommend using `Homebrew <http://brew.sh/>` to install a friendlier version of Python that doesn't require root permissions; see http://docs.python-guide.org/en/latest/starting/install/osx/ for instructions.

**Python version**
The toolkit is not yet compatible with Python 3. For now, make sure you're using Python 2.7 to install and run everything.


Notebook extensions
^^^^^^^^^^^^^^^^^^^
MDT will automatically install and enable the ``nbmolviz`` and ``widgetsnbextensions`` for Jupyter if they're not already enabled (these extensions provide the interactive molecular visualization framework). You can list the installed extensions by running

.. code-block:: bash

    jupyter nbextension list


They can be turned on, if necessary, by running:

.. code-block:: bash

    jupyter nbextension enable --python nbwidgetsextension
    jupyter nbextension enable --python nbmolviz


Optional python installs
========================
MDT relies on several open source python packages to provide various functions.
It's not necessary to have any of them installed on your own machine, but it can be convenient to
have them locally, especially for development

OpenBabel
^^^^^^^^^
MDT makes heavy use of `OpenBabel's <https://openbabel.org/>`_ file i/o routines. While powerful,
OpenBabel and its Python bindings can be challenging to compile. It's highly recommended to use
a package manager:

 * For Debian and Ubuntu, use ``apt-get install openbabel python-openbabel``
 * For MacOS X, use ``brew install --with-python open-babel``.

See `the OpenBabel docs <https://openbabel.org/docs/dev/Installation/install.html>`_ for more
information.

OpenMM
^^^^^^
We recommend the `pre-compiled binary installers available directly from the OpenMM
developers. <https://simtk.org/frs/?group_id=161>`_

PySCF
^^^^^
MDT uses quantum chemical implementations and basis set logic from the excellent `PySCF library <http://sunqm.net/pyscf/>`_. You can install it locally by running

.. code-block:: bash

    pip install git+https://github.com/sunqm/pyscf

If you run into problems, see the `documentation <http://sunqm.net/pyscf/>`_ and
`GitHub README. <https://github.com/sunqm/pyscf>`_.


Changing where your jobs run
============================

The toolkit is built to run jobs using the `Docker containerization technology <https://www.docker.com/>`_
(which *has nothing to do with molecular docking*).  Docker eliminates the need to configure or compile
software on different computers.

By default, MDT is configured to use a free cloud-based docker cluster provided by Autodesk
Research. If you'd like to run jobs on your local machine, you'll need to install a couple more
things.



Running locally with Docker
^^^^^^^^^^^^^^^^^^^^^^^^^^^
First, create or edit a file at ``$HOME/.moldesign/moldesign.yml`` with the line

.. code-block:: yaml

    engine_type: docker

Next, install Docker if necessary (version 1.11 or higher is required):

- *Mac*: Download and install `Docker for Mac <https://docs.docker.com/docker-for-mac/>`_.
- *Windows*: Download and install `Docker for Windows <https://docs.docker.com/docker-for-windows/>`_.
- *Linux*: `Follow the instructions for your distribution <https://docs.docker.com/engine/installation/linux/>`_.

Once Docker is up and running, make sure to allocate enough RAM - 4 GB will work well for the
included example jobs.

Running locally with CloudComputeCannon
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Our group has also developed
`CloudComputeCannon <https://www.npmjs.com/package/cloud-compute-cannon>`_, a lightweight,
Docker-based job scheduling system which is more suitable for production than a bare Docker engine.

You'll need Docker installed locally (see steps above). To install CCC:

1. `Install the Node.js <https://nodejs.org/en/>`_ javascript interpreter if necessary.
2. Update NPM if necessary: ``npm install npm -g``
3. Do a global install of cloud compute cannon: ``npm install -g cloud-compute-cannon``

To run it:

- To **start** the CCC scheduler, make sure Docker is running locally, then run ``ccc server-install``
- To **stop** the CCC scheduler, run ``ccc server-stop``

Finally, update your MDT configuration to point to the CCC server by default by putting these lines in
``$HOME/.moldesign/moldesign.yml``:

.. code-block:: yaml

    engine_type: cloudcomputecannon
    default_ccc_server: localhost:9000





System-specific installation
============================

Use these instructions to get the dependencies installed on your machine, then proceed with the
`basic installation <#basic-installation>`_.

Mac OS X
^^^^^^^^
Install `homebrew <http://brew.sh>`_ if it's not already installed:

.. code-block:: bash

   /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

Install Python:

.. code-block:: bash

   $ brew install python

Ubuntu / Debian
^^^^^^^^^^^^^^^
.. code-block:: bash

   apt-get install python python-pip

Red Hat / Fedora / CentOS
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   yum install python python-pip

Windows
^^^^^^^
Windows is not supported at this time, although if the Python dependencies can be installed, MDT should work fine. We highly recommend running your notebooks in Chrome or Firefox. We're working on Windows support - if you have any problems, please `open an issue on GitHub <https://github.com/autodesk/mol
design/issues>`_.