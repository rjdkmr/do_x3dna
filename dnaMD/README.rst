.. |3DNA package| raw:: html

   <a href="http://x3dna.org" target="_blank">3DNA package</a>

.. |dnaMD.py| raw:: html

   <a href="https://github.com/rjdkmr/do_x3dna/blob/master/Python_API/dnaMD.py" target="_blank">dnaMD.py</a>

.. |mitochondiral DNA| raw:: html

   <a href="http://pdb.org/pdb/explore/explore.do?structureId=3tmm" target="_blank">mitochondiral DNA</a>

.. |nucleosome| raw:: html

   <a href="http://pdb.org/pdb/explore/explore.do?structureId=3Ut9" target="_blank">nucleosome</a>

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

From the MD trajectory, several `output files <http://rjdkmr.github.io/do_x3dna/usage.html#output-files>`_ are generated.
For easy analysis of the obtained data, do_x3dna package includes a Python code |dnaMD.py|, which contains
several `methods <http://rjdkmr.github.io/do_x3dna/apidoc.html>`_ for the analysis of structural descriptors.

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

**For detailed documentation about the do_x3dna, please visit  `do_x3dna home-page<http://rjdkmr.github.io/do_x3dna>`_.**


Citations
---------

**Please cite the follwoing publications:**

* | Xiang-Jun Lu & Wilma K. Olson (2003)
  | 3DNA: a software package for the analysis, rebuilding and visualization of three-dimensional nucleic acid structures.
  | *Nucleic Acids Res.* 31(17), 5108-21.

* | Rajendra Kumar and Helmut Grubmüller (2015)
  | `do_x3dna: a tool to analyze structural fluctuations of dsDNA or dsRNA from molecular dynamics simulations <https://doi.org/10.1093/bioinformatics/btv190>`_
  | *Bioinformatics* (2015) 31 (15): 2583-2585.
