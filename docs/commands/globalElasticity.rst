globalElasticity
================

Global Elastic Properties of the DNA

This can be used to Global Elastic Properties of the DNA from the simulations.

**Usage:**

.. code-block:: bash

    usage: dnaMD globalElasticity [-h] [-i parameter.h5] [-o output.csv]
                                  [-ot elasicity_vs_time.csv]
                                  [-tbp total-bp-number] [-estype esType]
                                  [-bs bp/s-start-number] [-be bp/s-end-number]
                                  [-fgap frames-to-skip] [-paxis X] [-em block]
                                  [-gt gmx analyze] [-fbp 1]


**Optional arguments:**

.. code-block:: bash

    -h, --help            show this help message and exit
    -i parameter.h5, --input parameter.h5
                        Name of input file (hdf5 file).
                        This file should contain the required parameters. It should be hdf5 storage file.

    -o output.csv, --output output.csv
                        Name of output file with csv extension.
                        This file will contain the elasticity modulus matrix where values will be separated
                        by comma. Since modulus matrix is also shown as screen output, this option is not
                        necessary.

    -ot elasicity_vs_time.csv, --output-time elasicity_vs_time.csv
                        Calculate elasticity as a function of time and save in this csv format file.
                        It can be used to obtained elastic moduli as a function of time to check their
                        convergence. If this option is used, -fgap/--frame-gap is an essential option.

                        NOTE: Elastic properties cannot be calculated using a single frame because
                        fluctuations are required. Therefore, here time means trajectory between zero
                        time to given time.

    -tbp total-bp-number, --total-bp total-bp-number
                        Total number of basepair in DNA/RNA.
                        It is an essential input.

    -estype esType, --elasticity-type esType
                        Elastic Properties type.
                        Two keywords are accepted: "BST" or "ST".
                        * "BST" : Bending-Stretching-Twisting --- All motions are considered
                        * "ST"  : Stretching-Twisting --- Bending motions are ignored.

                        WARNING: For accurate calculation of bending motions, DNA structures in trajectory must
                        be superimposed on a reference structure (See Publication's Method Section).

    -bs bp/s-start-number, --bp-start bp/s-start-number
                        First BP/BPS of DNA after which parameter will be extracted.
                        If it is not given, first basepair or base-step will be considered.

    -be bp/s-end-number, --bp-end bp/s-end-number
                        Last BP/BPS of DNA upto which parameter will be extracted.

                        If it is not given, last basepair or base-step will be considered.

    -fgap frames-to-skip, --frame-gap frames-to-skip
                        Number of frames to skip for next calculation
                        When calculating elastic modulus as a function of time, this option will determine
                        the time-gap between each calculation point. Lower the time-gap, slower will be the
                        calculation.
    -paxis X, --principal-axis X
                        Principal axis parallel to global helical-axis
                        Three keywords are accepted: "X", "Y" and "Z". Only require when bending motions are
                        included in the calculation.

    -em block, --error-method block
                         Error of elastic modulus
                        If this option is used, elastic modulus will be calculated as a function of time. Therefore,
                        options such as frameGap will be essential.

                        Error methods are as following:
                        * "none" : No error calculation (Default).
                        * "acf": Using autocorrelation function to determine autocoprrelation time and used as time
                                 to get the independent frame.
                        * "block": Block averaging error
                        * "std": standard deviation

                        In case of "acf" and "block", gromacs tool "g_analyze" or "gmx analyze" will be used. Either
                        of these tools should be in path for error calculation.

    -gt gmx analyze, --gromacs-tool gmx analyze
                        Tools to calculate autocorrelation time or bloack averaging error.
                        By default it is g_analyze (Gromacs-4.5.x/4.6.x versions). For newer versions, use "gmx analyze".

    -fbp 1, --first-bp 1  Basepair number of first base-pair.
                        Usually it is one. Therefore, if this option is not provided, base-pair
                        numbering will start from one.

                        In rare cases, base-pair numbering might start with other number. In those
                        cases, use this option to start numbering of basepair from other number than
                        one.


Example
-------

1. See example `here <../global_elasticity.html#Bending-Stretching-Twist-Modulus>`_.
2. See next example `here <../global_elasticity.html#Stretching-Twist-Modulus>`_.
3. See next example `here <../global_elasticity.html#Convergence-in-Modulus>`_.