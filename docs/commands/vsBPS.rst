vsBPS
=====

**Average parameters as a function of base-pair/step**

This can be used to calculate average and error of the given parameter of either a specific base-pair/step
or over a DNA segment during the simulations.

**Usage:**

.. code-block:: bash

    usage: dnaMD vsBPS [-h] [-i L-BP_cdna.dat] [-o output.dat]
                       [-tbp total-bp-number] [-p parameter] [-em block]
                       [-mb bp/s number] [-gt g_analyze] [-bs bp/s-start-number]
                       [-be bp/s-end-number] [-mm sum-or-mean] [-fbp 1]


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

      -p parameter, --parameter parameter
                            Parameter name.
                            This parameter will be extracted from file. Ensure that parameter is present
                            in the file, otherwise wrong values will be extracted from file.

      -em block, --error-method block
                             Error method
                            It can be as following:
                            * "acf": Using autocorrelation function to determine autocoprrelation time and used as time
                                     to get the independent frame.
                            * "block": Block averaging error
                            * "std": standard deviation

                            In case of "acf" and "block", gromacs tool "g_analyze" or "gmx analyze" will be used. Either
                            of these tools should be in path for error calculation.

      -mb bp/s number, --merge-bps bp/s number
                            Number of consecutive base-pairs/steps to merge for creating the small non-overlapping DNA
                            segments. By default, averages and errors are calculated for each base-pair/step separately.
                            However, averages  and errors can also  be calculated for small non-overlapping DNA segment
                            by merging the parameters of consecutive base-pairs/steps.

      -gt g_analyze, --gromacs-tool g_analyze
                            Tools to calculate autocorrelation time or bloack averaging error.
                            By default it is g_analyze (Gromacs-4.5.x/4.6.x versions). For newer versions, use "gmx analyze".

      -bs bp/s-start-number, --bp-start bp/s-start-number
                            First BP/BPS of DNA after which parameter will be extracted.
                            If it is not given, first basepair or base-step will be considered.

      -be bp/s-end-number, --bp-end bp/s-end-number
                            Last BP/BPS of DNA upto which parameter will be extracted.
                            If it is not given, last basepair or base-step will be considered.

      -mm sum-or-mean, --merge-method sum-or-mean
                            Method to merge the parameters of consecutive basepairs/steps given by -mb/--merge-bps.

                            Currently accepted keywords are as follows:
                                * mean : Average of local parameters
                                * sum : Sum of local parameters

      -fbp 1, --first-bp 1  Basepair number of first base-pair.
                            Usually it is one. Therefore, if this option is not provided, base-pair
                            numbering will start from one.

                            In rare cases, base-pair numbering might start with other number. In those
                            cases, use this option to start numbering of basepair from other number than
                            one.


Example
-------

.. code-block:: bash

    dnaMD vsBPS -i fdna.h5 -o vsBPS.dat -tbp 60 -bs 5 -be 56 -p rise -mm sum -mb 4 -gt "gmx analyze" -em block

Following output is obtained in ``vsBPS.dat`` file.

::

    # bp(mid) 	 rise-avg 	 rise-error
    7 		     13.313 	 0.0161982
    11 		     13.5327 	 0.0174155
    15 		     13.5497 	 0.0191224
    19 		     13.376 	 0.0778458
    23 		     13.5107 	 0.0167237
    27 		     13.6694 	 0.0254001
    31 		     13.5624 	 0.0178427
    .
    .
    .

It can be plotted by xmgrace as following:

.. code-block:: bash

    xmgrace -settype xydy vsBPS.dat


The obtained plot is similar to the second plot
`shown here <../notebooks/base_steps_tutorial.html#Local-base-step-parameters-as-a-function-of-base-steps>`_
for bound DNA.
