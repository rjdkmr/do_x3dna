.. |3DNA package| raw:: html

   <a href="http://x3dna.org" target="_blank">3DNA package</a>


.. |mitochondiral DNA| raw:: html

   <a href="http://pdb.org/pdb/explore/explore.do?structureId=3tmm" target="_blank">mitochondiral DNA</a>


.. |do_x3dna homepage| raw:: html

  <a href="http://do-x3dna.readthedocs.io" target="_blank">do_x3dna homepage</a>


Introduction
============

do_x3dna has been developed for analysis of the DNA/RNA dynamics during the molecular dynamics simulations.
It consists of three main components:

    * `do_x3dna <http://do-x3dna.readthedocs.io/en/latest/do_x3dna_usage.html>`_
      --- To calculate structural descriptors of DNA/RNA from MD trajectory.

    * `dnaMD <http://do-x3dna.readthedocs.io/en/latest/dnaMD_usage.html>`_
      --- Command line tool to extract and analyze the data obtained from do_x3dna
      for **users without programming experiences**.

    * `dnaMD Python module <http://do-x3dna.readthedocs.io/en/latest/api_summary.html>`_
      --- To extract and analyze the data obtained from do_x3dna for
      **users with programming experiences**.


**Last Update: Aug. 2017**

**For detailed documentation about the do_x3dna, please visit:** |do_x3dna homepage|.

**For Questions and Discussions, please visit:** `do_x3dna forum <https://groups.google.com/forum/#!forum/do_x3dna>`_.

Release Note 2017
-----------------

  * `do_x3dna <do_x3dna_usage.html>`_ can be compiled and used with **GROMACS**
    **4.5.x**, **4.6.x**, **5.0.x**, **5.1.x**, and **2016.x** versions.

  * More user friendly --- `dnaMD <dnaMD_usage.html>`_ tools to analyze
    `do_x3dna data <do_x3dna_usage.html#output-files-table>`_ --- No
    programming experiences needed now to analyze do_x3dna data.

  * `Speed up dnaMD analysis with HDF5 file <using_hdf5.html>`_


Citations
---------

**Please cite the follwoing publications:**

* | Xiang-Jun Lu & Wilma K. Olson (2003)
  | 3DNA: a software package for the analysis, rebuilding and visualization of three-dimensional nucleic acid structures.
  | *Nucleic Acids Res.* 31(17), 5108-21.

* | Rajendra Kumar and Helmut Grubmüller (2015)
  | `do_x3dna: a tool to analyze structural fluctuations of dsDNA or dsRNA from molecular dynamics simulations <https://doi.org/10.1093/bioinformatics/btv190>`_
  | *Bioinformatics* (2015) 31 (15): 2583-2585.
