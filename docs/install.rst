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

Notebook extensions
^^^^^^^^^^^^^^^^^^^
MDT will automatically install and enable the ``nbmolviz`` and ``widgetsnbextensions`` for Jupyter if they're not already enabled (these extensions provide the interactive molecular visualization framework). You can list the installed extensions by running

.. code-block:: bash

    jupyter nbextension list


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

The toolkit is built to run jobs using the `Docker containerization technology <https://www.docker.com/>`_ (which *has nothing to do with molecular docking*).  Docker eliminates the need to configure or compile
software on different computers.

By default, MDT is configured to use a free cloud-based docker cluster provided by Autodesk
Research. If you'd like to run jobs on your local machine, you'll need to install a couple more
things.


Running jobs locally
--------------------

Using a docker-machine
^^^^^^^^^^^^^^^^^^^^^^
A recent version of Docker (>1.11) is required.

*Mac or Windows*: Download and install the `Docker Toolbox <https://www.docker
.com/products/docker-toolbox>`_.

*Linux*: `Follow the instructions for your distribution <https://docs.docker
.com/engine/installation/linux/>`_.

Next, create a docker-machine (ideally, it should have at least 4 GB of RAM and 40 GB of disk
space):

.. code-block:: bash

    $ docker-machine create --driver virtualbox --virtualbox-memory "4096" --virtualbox-disk-size "40000"


Running jobs on AWS
--------------------
coming soon

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