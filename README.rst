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

    * do_x3dna
    * dnaMD
    * dnaMD Python module

do_x3dna
--------
It is the wrapper tool, which uses |3DNA package| to calculate several structural
descriptors of DNA/RNA from the GROMACS MD trajectory. It executes 3DNA tools to
calculate these descriptors and subsequently, extracts these output and saves in to
external output files as a function of time.

From the MD trajectory, several `output files <http://do-x3dna.readthedocs.io/en/latest/do_x3dna_usage.html#output-files>`_ are generated.
For easy analysis of the obtained data, do_x3dna package includes a `dnaMD`_ tool and Python module, which contains
several `tools <http://do-x3dna.readthedocs.io/en/latest/dnaMD_usage.html>`_ and
 `methods <http://do-x3dna.readthedocs.io/en/latest/api_summary.html>`_ for the analysis of structural descriptors.

**GROMACS versions supported:** 4.5.x, 4.6.x, 5.0.x, 5.1.x, and 2016.x versions.
For higher versions, a PDB file can be used in place of a GROMACS **tpr** file.

.. note::
    do_x3dna can be used with trajectory files that are obtained from other MD packages such as NAMD and AMBER.
    Input trajectory files should be converted in to Gromacs format trajectory files. A PDB file could be used in place
    of a GROMACS **tpr** file.

dnaMD
-----
`do_x3dna`_ generates large amount of data and analyzing these data is difficult.
Therefore, ``dnaMD`` is developed to analyze data obtained from `do_x3dna`_. It contains
set of tools to analyze the data.


dnaMD Python module
-------------------
`dnaMD`_ is written in Python and it is also available as Python module.
It can be used in Python scripting for fast and flexible analysis of `do_x3dna`_
data.

**Last Update: August. 2017**

**For detailed documentation about the do_x3dna, please visit |do_x3dna homepage| .**


Citations
---------

**Please cite the follwoing publications:**

* | Xiang-Jun Lu & Wilma K. Olson (2003)
  | 3DNA: a software package for the analysis, rebuilding and visualization of three-dimensional nucleic acid structures.
  | *Nucleic Acids Res.* 31(17), 5108-21.

* | Rajendra Kumar and Helmut Grubmüller (2015)
  | `do_x3dna: a tool to analyze structural fluctuations of dsDNA or dsRNA from molecular dynamics simulations <https://doi.org/10.1093/bioinformatics/btv190>`_
  | *Bioinformatics* (2015) 31 (15): 2583-2585.
