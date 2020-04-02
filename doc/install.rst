.. _installation:

############
Installation
############

Download
========

Go to the :doc:`download` page. Download LFA Lab from this page and unpack it
when necessary.

Prepare your System
===================

LFA Lab depends on the following packages:

- C++ compiler
- `CMake <http://www.cmake.org/>`_
- `Eigen3 <http://eigen.tuxfamily.org/>`_
- `Python3 <http://www.python.org/>`_
- `Swig <http://swig.org/>`_
- `NumPy <http://www.numpy.org/>`_
- `six <https://pypi.org/project/six/>`_
- Only for pre C++11 compilers: `Boost <http://www.boost.org/>`_

Furthermore, the following packages provide extra functionality, but are
*optional*:

- `Doxygen <http://www.doxygen.org/>`_
- `Sphinx <http://www.sphinx-doc.org/>`_
- `LAPACK <http://www.netlib.org/lapack/>`_
- `ARPACK <http://github.com/opencollab/arpack-ng/>`_
- `googletest <http://code.google.com/p/googletest/>`_
- `PIP <https://pip.pypa.io/en/stable/>`_

You can find the commands to install the dependencies for different operating
systems below.

Debian
------

To install the dependencies on Debian, run the following command::

  sudo apt-get install -y \
      g++ cmake python3-dev python3-numpy python-six swig \
      libeigen3-dev liblapack-dev python-matplotlib

Fedora
------

To install the dependencies on Fedora, run the following command::

  sudo dnf install -y \
      gcc-c++ cmake python3-devel python3-numpy python-six swig \
      eigen3-devel lapack-devel python3-matplotlib

Mac OS X
--------

First, you need to install a current version of Python. Visit the
`Python Website <https://www.python.org/>`_ and download a recent version
of Python 3, which you then install. Make sure that this Python version
is your path by running the `Update Shell Profile.command`. After running
this command, you need to close all terminal windows.

Next, we need to install the C++ dependencies, needed to compile the native
Python extension of LFA Lab. An easy way to install these dependencies is
the `Homebrew package manager <http://brew.sh>`_ for Mac. Make sure you have
Homebrew installed. Then, change into the source directory of LFA Lab and
execute the following command::

  brew bundle

Then, you can use `PIP`_ to install the remaining (Python) dependencies::

  pip3 install --user --upgrade -r requirements.txt

*Warning*: You might have multiple Python versions on your machine. You have
to make sure that all dependencies and LFA Lab are installed with the same
version.

Quick Install
=============

After you have prepared your system, you can build and install LFA Lab using
pip. Change into the source directory of LFA Lab and execute::

  pip3 -mpip install --user .

To check your installation, open a *new* terminal window and run::

  python3 -mlfa_lab

You should see a message, which tells you that LFA Lab has been installed
successfully.

If you want to customize the installation of in case the build fails, see the
manual installation below.

.. _build_lfa_lab:

Build LFA Lab
=============

Execute the following commands in the source directory of LFA Lab::

    cmake [OPTIONS] .
    make

(Do not forget the dot at the end.)

In case the build fails or you want to tweak your installation you can use the
following options.

Options
-------

Install into a per-user directory::

    -DUSER_INSTALL=ON

Turns the compiler optimization on::

    -DCMAKE_BUILD_TYPE=Release

Choose to use LAPACK for certain operations. If a good LAPACK/BLAS
implementation is available, this will speed up the program essentially::

    -DWITH_LAPACK=[ON|OFF]

Choose to use ARPACK. This will speed up the program if large spectra
need to be analyzed. Arpack, however, might not be able to compute the spectra
for certain tricky problems::

    -DWITH_ARPACK=[ON|OFF]

Set other prefices wich will be searched. For example if you installed
some of the libraries in $HOME/.local run::

    -DCMAKE_PREFIX_PATH=/other/prefix1;/other/prefix2

For example::

    cmake -DCMAKE_PREFIX_PATH=$HOME/.local [OTHER OPTIONS] .

Documentation
-------------

To build the documentation you can run::

    make sphinx-doc

This command requires `Sphinx`_.

The C++-Core modules can be documented using::

    make doxygen

Installation
============

To install LFA Lab just run::

    sudo make install

If you just want to use the software without installation, you can run::

    source setup-env.py

instead. This command will setup the current shell session such that you can
use LFA Lab.

You can now use LFA Lab. Take a look at the :doc:`tutorial` page to find out
how to use it.

