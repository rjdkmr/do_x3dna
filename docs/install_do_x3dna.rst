Download and Installation of do_x3dna for GROMACS
=================================================

Download
--------

**do_x3dna package can be downloaded by either of following two steps.**

1. Use ``git``:

::

    git clone -b master --single-branch https://github.com/rjdkmr/do_x3dna


2. Download either `zipball <https://github.com/rjdkmr/do_x3dna/archive/master.zip>`_ or `tarball <https://github.com/rjdkmr/do_x3dna/archive/master.tar.gz>`_.

Installation
------------

do_x3dna requires one of the following versions of Gromacs:

    * Gromacs-4.5.x
    * Gromacs-4.6.x
    * Gromacs-5.0.x
    * Gromacs-5.1.x
    * Gromacs-2016.x
    * Gromacs-2018.x

For Gromacs **4.5.x** and **4.6.x** versions libraries such as ``libmd.a/.so``,
``libgmx.a/.so`` and ``libgmxana.a/.so`` are required.

For Gromacs **5.0.x**, **5.1.x**, **2016.x** and **2018.x** version ``libgromacs.a/.so`` is
required.


Steps
-----
Follow these steps to install do_x3dna

.. code-block:: bash

    cd do_x3dna
    mkdir build
    cd build
    cmake ..  -DGMX_PATH=/opt/gromacs -DCMAKE_INSTALL_PREFIX=/opt/do_x3dna
    make
    sudo make install


Directory ``/opt/gromacs`` should contains ``include`` and ``lib`` directories.
Gromacs version will be automatically detected and accordingly processed during
compilations.

  * ``-DGMX_PATH`` : path to GROMACS installation directory.

  * ``-DCMAKE_INSTALL_PREFIX``: path to a directory where do_x3dna will be installed.
    If this option is not provided, do_x3dna will install either in ``/usr/local/bin``
    or in ``$HOME/bin`` directory.

  * Instead of ``GMX_PATH``, ``CMAKE_PREFIX_PATH`` can be used before running cmake as follows:

      * ``export CMAKE_PREFIX_PATH=/opt/gromacs``
      * ``export CMAKE_PREFIX_PATH=/opt/gromacs:/path/to/fftw``

  * In case of **ERROR**: **FFTW library file libfftw3f.so or libfftw3f.a not found...**
    use ``-DFFTW_LIB=/path/to/fftw/lib/`` to give a path for the directory containing either of the two files.
