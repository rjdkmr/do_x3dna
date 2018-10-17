localEnergy
===========

Local Deformation Energy of the DNA

This can be used to Local Deformation Energy of the DNA from the simulations.
At first, elastic matrix from reference DNA (most often free or unbound DNA)
is calculated and subsequently this matrix is used to calculate local deformation
energy of probe DNA (most often bound DNA).

WARNING: The option "-o/--output" cannot be used with "-os/--output-segments"
         simultaneously due to incompatible usage of "-bs/--bps-start" and
         "-be/--bps-end".

WARNING: In case of "-o/--output" difference between "-bs/--bps-start"
         and "-be/--bps-end" should be less than 4 while "--span" should be
         less than 4 in case of "-os/--output-segments".


**Usage:**

.. code-block:: bash

    usage: dnaMD localEnergy [-h] [-ir ref_dna.h5] [-ip probe_dna.h5]
                             [-o energy_vs_time.csv] [-os segments_elasticity.csv]
                             [-eu kJ/mol] [-et "full,diag,slide,twist"]
                             [-tbp total-bp-number] [-helical]
                             [-bs basepair-step-start-number]
                             [-be basepair-step-end-number] [-span segment-span]
                             [-em block] [-gt gmx analyze] [-fbp 1]

**Optional arguments:**

.. code-block:: bash

      -h, --help            show this help message and exit
      -ir ref_dna.h5, --input-ref ref_dna.h5
                            Name of input reference file (hdf5 file).
                            File containing parameters of reference DNA for which local elastic properties will
                            be calculated. Most often it is free or unbound DNA.

                            This file should contain the required parameters. It should be hdf5 storage file.

      -ip probe_dna.h5, --input-probe probe_dna.h5
                            Name of input probe file (hdf5 file).
                            File containing parameters of probe DNA for which local deformation energy will
                            be calculated. Most often it is bound DNA.

                            This file should contain the required parameters. It should be hdf5 storage file.

      -o energy_vs_time.csv, --output energy_vs_time.csv
                            Calculate energy as a function of time and save in this csv format file.
                            This file will contain the energy values of a local small DNA segment as a
                            function of time.

      -os segments_elasticity.csv, --output-segments segments_elasticity.csv
                            Calculate energy of consecutive overlapped DNA segments and save in this csv format file.
                            It enables the calculation of local deformation energies of small overlapped DNA segments with
                            error. This file will contain the energies of small overlapped DNA segments of
                            length given by "-s/--span". Error will be also calculated along with the average values.

      -eu kJ/mol, --energy-unit kJ/mol
                            Unit of energy. (Default: kJ/mol)

                            Allowed units are:
                             * kT
                             * kJ/mol
                             * kcal/mol

      -et "full,diag,slide,twist", --energy-terms "full,diag,slide,twist"
                            Energy terms to be calculated.
                            For which motions, energy should be calculated.

                            Following keywords are available:
                                * 'all'                   : (Default) All below listed energy terms will be calculated
                                * 'full'                  : Use entire elastic matrix -- all parameters with their coupling
                                * 'diag'                  : Use diagonal of elastic matrix -- all motions but no coupling
                                * 'shift' or 'x-disp'
                                * 'slide' or 'y-idsp'
                                * 'rise'  or 'h-rise'
                                * 'tilt'  or 'inclination'
                                * 'roll'  or 'tip'
                                * 'twist' or 'h-twist'

                            The terms should provided as comma separated values. e.g. -et "full,diag,shift,slide,twist".

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

                            In case of "-o/--output" this is firs bp-s of the segment while for
                            "-os/--output-segments", it is first bp-s of first overlapped segments.
      -be basepair-step-end-number, --bps-end basepair-step-end-number
                            Last basepair-step (bp-s) of DNA upto which parameter will be extracted.

                            In case of "-o/--output" this is last bp-s of the segment while for
                            "-os/--output-segments", it is last bp-s of last overlapped segments.

      -span segment-span, --span segment-span
                            Length of overlapping (local) DNA segments.
                            It is essential when "-os/--output-segments" is used. It should not be larger than four.

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

1. See example `here <../global_elasticity.html#Local-deformation-energy-of-a-local-small-segment>`_.
2. See next example `here <../global_elasticity.html#Deformation-energy-of-the-consecutive-overlapped-DNA-segments>`_.