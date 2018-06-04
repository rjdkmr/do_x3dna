histogram
=========

**Parameter distribution during simulation**

This can be used to calculate distribution of a parameter of either a specific base-pair/step
or over a DNA segment during the simulations.

**Usage:**

.. code-block:: bash

    usage: dnaMD histogram [-h] [-i L-BP_cdna.dat] [-o output.dat]
                           [-tbp total-bp-number] [-bins bins] [-p parameter]
                           [-bs bp/s-start-number] [-be bp/s-end-number]
                           [-mm sum-or-mean] [-fbp 1]


**Optional arguments:**

.. code-block:: bash

    -h, --help            show this help message and exit
    -i L-BP_cdna.dat, --input L-BP_cdna.dat
                        Name of input file (from do_x3dna or hdf5 file).
                        This file should contain the required parameters. It can be a file either
                        produced from do_x3dna or hdf5 storage file.

    -o output.dat, --output output.dat
                        Name of output file.
                        The extracted output will be written in output file.

    -tbp total-bp-number, --total-bp total-bp-number
                        Total number of basepair in DNA/RNA.
                        It is an essential input.

    -bins bins, --bins bins
                         Number of bins in the histogram
                        Default value is 30.

    -p parameter, --parameter parameter
                        Parameter name.
                        This parameter will be extracted from file. Ensure that parameter is present
                        in the file, otherwise wrong values will be extracted from file.

    -bs bp/s-start-number, --bp-start bp/s-start-number
                        First BP/BPS of DNA after which parameter will be extracted.
                        If it is not given, first basepair or base-step will be considered.

    -be bp/s-end-number, --bp-end bp/s-end-number
                        Last BP/BPS of DNA upto which parameter will be extracted.

                        If it is not given, parameter for only a single bp/s given with -bs/--bp-start
                        option will be extracted.

    -mm sum-or-mean, --merge-method sum-or-mean
                        Method to merge the parameter of a DNA segment from local parameters
                        of all base-pairs/steps that are within the range given by '-bs' and '-be'.

                        Currently accepted keywords are as follows:
                            * mean : Average of local parameters
                            * sum : Sum of local parameters

                        When only "-bs" option is provided without "-be", then -mm/--merge-method is
                        not required.

    -fbp 1, --first-bp 1  Basepair number of first base-pair.
                        Usually it is one. Therefore, if this option is not provided, base-pair
                        numbering will start from one.

                        In rare cases, base-pair numbering might start with other number. In those
                        cases, use this option to start numbering of basepair from other number than
                        one.


Example
-------

.. code-block:: bash

    dnaMD histogram -i pdna.h5 -o histogram.dat -tbp 60 -bs 20 -be 45 -p rise -mm sum -bins 20

Following output is obtained in ``histogram.dat`` file.

::

    # "rise" 	 Density
    84.082	0.0
    84.278	0.00509694
    84.474	0.0
    84.67	0.00509694
    84.866	0.00509694
    85.062	0.0203878
    85.258	0.0356786
    85.454	0.107036
    85.65	0.0560664
    .
    .
    .

It can be plotted by xmgrace as following:

.. code-block:: bash

    xmgrace histogram.dat


The obtained plot is similar to the histogram plot
`shown here <../notebooks/base_steps_tutorial.html#Distribution-of-local-base-steps-parameters-during-MD-simulations>`_
for 20-45 bp length bound DNA.
