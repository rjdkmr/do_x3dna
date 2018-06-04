localDeformation
================

**Deformation of local parameters in probe DNA with respect to a reference DNA**

This can be used to calculate deformation in the given parameters of a probe DNA with respect to a reference DNA along
the base-pairs/steps.

.. note:: Number of input segments/bp/bps should match between probe and reference DNA. It means ([-bePrb] - [-bsPrb])
          should be equal to ([-beRef] - [-bsRef]).

**Usage:**

.. code-block:: bash

    usage: dnaMD localDeformation [-h] [-ir L-BP_fdna.dat] [-ip L-BP_cdna.dat]
                                  [-o output.dat] [-oxb ref]
                                  [-tbpPrb total-bp-number]
                                  [-tbpRef total-bp-number] [-p parameter]
                                  [-em block] [-mb bp/s number] [-gt g_analyze]
                                  [-bsPrb bp/s-start-number]
                                  [-bePrb bp/s-end-number]
                                  [-bsRef bp/s-start-number]
                                  [-beRef bp/s-end-number] [-mm sum-or-mean]
                                  [-fbpPrb 1] [-fbpRef 1]


**Optional arguments:**

.. code-block:: bash

    -h, --help            show this help message and exit
    -ir L-BP_fdna.dat, --input-reference L-BP_fdna.dat
                        Name of input file for reference DNA (from do_x3dna or hdf5 file).
                        This file should contain the required parameters. It can be a file either
                        produced from do_x3dna or hdf5 storage file.

    -ip L-BP_cdna.dat, --input-probe L-BP_cdna.dat
                        Name of input file for probe DNA (from do_x3dna or hdf5 file).
                        This file should contain the required parameters. It can be a file either
                        produced from do_x3dna or hdf5 storage file.

    -o output.dat, --output output.dat
                        Name of output file.
                        The extracted output will be written in output file.

    -oxb ref, --xaxis-bp-source ref
                         Basepair/step number on the x-axis in output file
                        To write the base-pair/step number of either probe or reference DNA in the output file.

                        If base-pair numbering is different in reference and probe DNAs, then this option can be used to
                        write the base-pair number of either probe (-oxb probe) or reference (-oxb ref) DNA in the output file.

    -tbpPrb total-bp-number, --total-bp-probe total-bp-number
                        Total number of basepair in probe DNA.
                        It is an essential input.

    -tbpRef total-bp-number, --total-bp-Ref total-bp-number
                        Total number of basepair in reference DNA.
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

    -bsPrb bp/s-start-number, --bp-start-probe bp/s-start-number
                        First BP/BPS of probe DNA after which parameter will be extracted.
                        If it is not given, first basepair or base-step will be considered.

    -bePrb bp/s-end-number, --bp-end-probe bp/s-end-number
                        Last BP/BPS of probe DNA upto which parameter will be extracted.
                        If it is not given, last basepair or base-step will be considered.

    -bsRef bp/s-start-number, --bp-start-ref bp/s-start-number
                        First BP/BPS of reference DNA after which parameter will be extracted.
                        If it is not given, first basepair or base-step will be considered.

    -beRef bp/s-end-number, --bp-end-ref bp/s-end-number
                        Last BP/BPS of reference DNA upto which parameter will be extracted.
                        If it is not given, last basepair or base-step will be considered.

    -mm sum-or-mean, --merge-method sum-or-mean
                        Method to merge the parameters of consecutive basepairs/steps given by -mb/--merge-bps.

                        Currently accepted keywords are as follows:
                            * mean : Average of local parameters
                            * sum : Sum of local parameters

    -fbpPrb 1, --first-bp-probe 1
                        Basepair number of first base-pair in probe DNA.
                        Usually it is one. Therefore, if this option is not provided, base-pair
                        numbering will start from one.

    -fbpRef 1, --first-bp-ref 1
                        Basepair number of first base-pair in reference DNA.
                        Usually it is one. Therefore, if this option is not provided, base-pair
                        numbering will start from one.



Example
-------

.. code-block:: bash

    dnaMD localDeformation -ir pdna.h5 -ip fdna.h5 -o deformation.dat -tbpPrb 60 -tbpRef 60 -bsPrb 5 -bePrb 56 \
    -bsRef 5 -beRef 56 -mm sum -mb 4 -em block -p shift -gt "gmx analyze"


Following output is obtained in ``deformation.dat`` file.

::

    # bp(mid) 	 shift-avg 	 shift-error
    7		0.0539361	0.0486333
    11		0.122707	0.0485817
    15		0.206374	0.269953
    19		-0.0630569	0.228433
    23		-0.535754	0.247015
    27		0.697832	0.253653
    31		-0.453067	0.0476092
    .
    .
    .

It can be plotted by xmgrace as following:

.. code-block:: bash

    xmgrace -settype xydy deformation.dat


The obtained plot is similar to the second plot
`shown here <../notebooks/base_steps_tutorial.html#Deviation-in-parameters-of-bound-DNA-with-respect-to-free-DNA>`_
for bound DNA.
