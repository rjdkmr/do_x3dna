do_x3dna for VMD
================
do_x3dna can be also used as VMD plugin or Tcl package with VMD. It is very similar to original do_x3dna tool,
which is executed with GROMACS format files.

.. note:: do_x3dna VMD plugin does **not depend** on GROMACS and original do_x3dna tool. It is standalone and independent plugin.

Following parameters are calculated using 3DNA package for each frame and written in separate files as a function of time.

* Base-pairs
* Hydrogen bonds between base-pairs
* Base-pairs parameters
* Base-steps parameters
* Helical Base-steps parameters
* Local helical axis
* Major and minor grooves
* Local helical radius
* Backbone dihedral angles (``alpha``, ``beta``, ``gamma``, ``delta``, ``epsilon``, ``zeta`` and ``chi``)
* Sugar dihedral angles (``v0``, ``v1``, ``v2``, ``v3`` and ``v4``)

These files could be used with the ``dnaMD`` module and tools for further analysis.

Download
--------

**do_x3dna package can be downloaded by either of following two steps.**

1. Use ``git``:

::

    git clone -b master --single-branch https://github.com/rjdkmr/do_x3dna


2. Download either `zipball <https://github.com/rjdkmr/do_x3dna/archive/master.zip>`_ or `tarball <https://github.com/rjdkmr/do_x3dna/archive/master.tar.gz>`_.


Installation
------------
Follow these steps to install do_x3dna as VMD plugin.

.. code-block:: bash

    cd do_x3dna
    cd vmd
    cp -r do_x3dna /usr/local/lib/vmd/plugins/noarch/tcl/.


If VMD is not at default location, copy ``vmd/do_x3dna`` in ``plugins/noarch/tcl/`` vmd directory.

How to use?
-----------

To use do_x3dna VMD plugin, 3DNA package should be installed and ``$X3DNA`` environment
variable (Detail is given in 3DNA manual) should be defined.

It can be either used from VMD TkConsole or through a Tcl scrip.

For example, below section can be used directly in VMD TkConsole or save it as a tcl script file.

.. code-block:: tcl

  package require do_x3dna

  set trajectories {traj-1.dcd traj-2.dcd}
  do_x3dna -tp example.psf -refFrame example.pdb -trajs $trajectories -dnaAtoms "nucleic" -fitAtoms "name P"

  # Uncomment it when using as script
  # quit


Above example can be run through a script.

.. code-block:: bash

  vmd -dispdev text -e run_script.tcl



.. note:: We recommend to use do_x3dna plugin through script with VMD as shown above.
          ``-dispdev text`` execute VMD in text mode as command and does not open
          Display window.

**Execute following command to get full help**

.. code-block:: tcl

    package require do_x3dna
    do_x3dna -help



Usage
-----

.. note:: Before running do_x3dna, make sure 3DNA package is correctly working for
          the input DNA/RNA structure. To check it, run do_x3dna with only reference frame
          instead of trajectory file.
          ``-noisy`` option switch-on the output message from the 3DNA package on the display
          that would be necessary to troubleshoot problems related to the 3DNA package
          Most common problem is the residue names mismatch in input DNA/RNA structure
          and the 3DNA package dictionary.

.. note:: Only PBC corrected trajectory files should be used as inputs.

.. note:: If ``-fitAtoms`` is provided, this atoms group is used to superimpose current
          frame on reference frame. Most of the parameters are unaffected by this fitting, however
          coordinates of the local helical axis could mismatch with the input coordinates with the
          DNA/RNA.

.. note:: By default, `bigdcd <http://www.ks.uiuc.edu/Research/vmd/script_library/scripts/bigdcd/>`_
          routine is used because it enables the trajectory reading frame by frame and allows to use
          very large trajectory. It can be switched-off by using ``-loadAll`` option.

Summary
-------

.. code-block:: tcl

  do_x3dna -help -tp <file> -ftp <file type>
           -dnaAtoms <DNA atoms> -fitAtoms <Atoms to fit>
           -refFrame <Reference file> -refFrameType <reference file type>
           -trajs <List of trajectory files> -trajType <trajectory file type>
           -suffix <output suffix> -dt <time-step ps>
           -noisy -loadAll

OPTIONS
-------

``-help``
~~~~~~~~~
**Optional**

Show help and exit.

``-tp``
~~~~~~~
**Required**

Topology-Parameter file with or without coordinates (psf/prmtop/pdb/gro).
If it is pdb/gro containing coordinates of atoms, this will be used as
reference frame for fitting and base-pair calculation. In case of coordinate
files such as pdb/gro, ``-refFrame`` will be only considered as an additional
frame as a part of trajectory.

``-ftp``
~~~~~~~~
**Optional**

File format of Topology-Parameter file. If not provided, file format
is automatically determined by VMD using extension.

``-dnaAtoms``
~~~~~~~~~~~~~
**Optional**

Atoms belonging to DNA/RNA. If not provided, by default ``"nucleic"`` is
used to select atoms.

``-fitAtoms``
~~~~~~~~~~~~~
**Optional**
Atoms for fitting to reference frame. If not provided, fitting will be
not done.

``-refFrame``
~~~~~~~~~~~~~
**Required/Optional**

Reference coordinate file for fitting and base-pair calculation. In case of psf
and prmtop file, this is required because these files do not contain
coordinates.

``-refFrameType``
~~~~~~~~~~~~~~~~~
**Optional**

File format of Reference frame file. If not provided, file format
is automatically determined by VMD using extension.

``-trajs``
~~~~~~~~~~
**Optional**

List of trajectory files. If not provided, only first frame will be used
for calculation. It is useful to check whether 3DNA is working correctly
with ``-noisy`` option.

``-trajType``
~~~~~~~~~~~~~
**Optional**

File format of trajectory file. If not provided, file format
is automatically determined by VMD using extension.

``-suffix``
~~~~~~~~~~~
**Optional**

Suffix for all output file names. If not provided, ``"g"`` is used by default.

``-dt``
~~~~~~~
**Optional**

Time-step between frames in ps. It is dumped in output files and make
easier for plotting. If not provided, by default, its value is 1 ps.

``-noisy``
~~~~~~~~~~
**Optional**

Display messages and errors from 3DNA programs ``find_pair`` and ``analyze``.

``-loadAll``
~~~~~~~~~~~~
**Optional**

With this option, at first all trajectory will be loaded and then 3DNA
will be executed for each frame. Therefore, large memory will be consumed
for large trajectory.


Output Files
------------

Following files are generated from ``do_x3dna``:

.. list-table:: List of output files from do_x3dna
    :widths: 1, 4
    :header-rows: 1
    :name: output-files-table-vmd

    * - File name
      - Output contents

    * - base_pairs_g.dat
      - Base-pairs

    * - h-bond_g.dat
      - Hydrogen bonds between base-pairs

    * - L-BP_g.dat
      - Base-pairs parameters

    * - L-BPS_g.dat
      - Base-steps parameters

    * - L-BPH_g.dat
      - Helical Base-steps parameters

    * - HelAxis_g.dat
      - Local helical axis coordinates

    * - MGroove_g.dat
      - Major and Minor grooves

    * - HelixRad_g.dat
      - Local helical radius

    * - BackBoneCHiDihedrals_g.dat
      - Backbone dihederal angles including Chi-dihedral

    * - SugarDihedrals_g.dat
      - Sugar dihedral angles including puckering type


Name of these files could be change by setting different suffix instead of ``g`` using ``-suffix`` option. These
files could be used with the ``dnaMD`` module and tools for further analysis.
