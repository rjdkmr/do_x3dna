localElasticity
===============
Local Elastic Properties of the DNA

This can be used to calculate local elastic properties of the DNA from the simulations. Here
local DNA segment referred to less than 5 base-pair (4 basepair-step)long.

WARNING: The option "-o/--output" and "-ot/--output-time" cannot be used with "-os/--output-segments"
         simultaneously due to incompatible usage of "-bs/--bps-start" and "-be/--bps-end". Both
         "-o/--output" and "-ot/--output-time" options can be used simultaneously.

WARNING: In case of "-o/--output" and "-ot/--output-time" difference between "-bs/--bps-start"
         and "-be/--bps-end" should be less than 4 while "--span" should be less than 4 in
         case of "-os/--output-segments".


**Usage:**

.. code-block:: bash

    usage: dnaMD localElasticity [-h] [-i parameter.h5] [-o output.csv]
                                 [-ot elasicity_vs_time.csv]
                                 [-os segments_elasticity.csv]
                                 [-tbp total-bp-number] [-helical]
                                 [-bs basepair-step-start-number]
                                 [-be basepair-step-end-number]
                                 [-span segment-span] [-fgap frames-to-skip]
                                 [-em block] [-gt gmx analyze] [-fbp 1]

**Optional arguments:**

.. code-block:: bash

      -h, --help            show this help message and exit
      -i parameter.h5, --input parameter.h5
                            Name of input file (hdf5 file).
                            This file should contain the required parameters. It should be hdf5 storage file.

      -o output.csv, --output output.csv
                            Name of output file with csv extension.
                            This file will contain the local elasticity matrix where values will be separated by comma.
                            This is also shown as screen output.

      -ot elasicity_vs_time.csv, --output-time elasicity_vs_time.csv
                            Calculate elasticity as a function of time and save in this csv format file.
                            It can be used to obtained local elasticity as a function of time to check their
                            convergence. If this option is used, -fgap/--frame-gap is an essential option.

                            NOTE: Elastic properties cannot be calculated using a single frame because
                            fluctuations are required. Therefore, here time means trajectory between zero
                            time to given time.

      -os segments_elasticity.csv, --output-segments segments_elasticity.csv
                            Calculate local elasticity of consecutive overlapped DNA segments and save in this csv format file.
                            It enables the calculation of local elastic properties of small overlapped DNA segments with
                            error. This file will contain the local elastic properties of small overlapped DNA segments of
                            length given by "-s/--span". As error will be calculated for each segment, '-fgap/--frame-gap'
                            is required.

      -tbp total-bp-number, --total-bp total-bp-number
                            Total number of basepair in DNA/RNA.
                            It is an essential input.

      -helical, --helical   Enable helical base-step parameters
                            By default, elasticity for base step-parameters (shift, slide, rise,
                            tilt, roll, twist) are calculated. If this option is used, helical
                            base-step parameters (x-disp, y-disp, h-rise, inclination, tip, h-twist)
                            will be used for elasticity.

      -bs basepair-step-start-number, --bps-start basepair-step-start-number
                            First basepair-step (bp-s) of DNA after which parameter will be extracted.
                            If it is not given, first bp-s will be considered.

                            In case of "-o/--output" and "-ot/--output-time" this is firs bp-s of
                            the segment while for "-os/--output-segments", it is first bp-s of
                            first overlapped segments.

      -be basepair-step-end-number, --bps-end basepair-step-end-number
                            Last basepair-step (bp-s) of DNA upto which parameter will be extracted.

                            In case of "-o/--output" and "-ot/--output-time" this is last bp-s of
                            the segment while for "-os/--output-segments", it is last bp-s of
                            last overlapped segments.

      -span segment-span, --span segment-span
                            Length of overlapping (local) DNA segments.
                            It is essential when "-os/--output-segments" is used. It should not be larger than four.

      -fgap frames-to-skip, --frame-gap frames-to-skip
                            Number of frames to skip for next calculation
                            When calculating elastic modulus as a function of time, this option will determine
                            the time-gap between each calculation point. Lower the time-gap, slower will be the
                            calculation.
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

1. See example `here <../global_elasticity.html#Local-elastic-properties>`_.
2. See next example `here <../global_elasticity.html#Local-elastic-properties-as-a-function-of-time>`_.
3. See next example `here <../global_elasticity.html#Local-elastic-properties-of-the-consecutive-overlapped-DNA-segments>`_.