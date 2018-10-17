globalEnergy
============

Global Deformation Energy of the DNA

This can be used to Global Deformation Energy of the DNA from the simulations. At first, elastic matrix from reference
DNA (most often free or unbound DNA) is calculated and subsequently this matrix is used to calculate deformation free
energy of probe DNA (most often bound DNA).


**Usage:**

.. code-block:: bash

    usage: dnaMD globalEnergy [-h] [-ir ref_dna.h5] [-ip probe_dna.h5]
                              [-o output.dat] [-et "full,diag,strecth,twist"]
                              [-tbp total-bp-number] [-estype esType]
                              [-bs bp/s-start-number] [-be bp/s-end-number]
                              [-paxis X] [-em block] [-gt gmx analyze] [-fbp 1]

**Optional arguments:**

.. code-block:: bash

      -h, --help            show this help message and exit
      -ir ref_dna.h5, --input-ref ref_dna.h5
                            Name of input reference file (hdf5 file).
                            File containing parameters of reference DNA for which global elastic properties will
                            be calculated. Most often it is free or unbound DNA.

                            This file should contain the required parameters. It should be hdf5 storage file.

      -ip probe_dna.h5, --input-probe probe_dna.h5
                            Name of input probe file (hdf5 file).
                            File containing parameters of probe DNA for which global deformation energy will
                            be calculated. Most often it is bound DNA.

                            This file should contain the required parameters. It should be hdf5 storage file.

      -o output.dat, --output output.dat
                            Name of output file in csv format.
                            This file will contain the energy values as a function of time.

      -et "full,diag,strecth,twist", --energy-terms "full,diag,strecth,twist"
                            Energy terms to be calculated.
                            For which motions, energy should be calculated.

                            Following keywords are available:
                            * all         : (Default) All below listed energy terms will be calculated
                            * full        : Use entire elastic matrix -- all motions with their coupling
                            * diag        : Use diagonal of elastic matrix -- all motions but no coupling
                            * b1          : Only bending-1 motion
                            * b2          : Only bending-2 motion
                            * stretch     : Only stretching motion
                            * twist       : Only Twisting motions
                            * st_coupling : Only stretch-twist coupling motion
                            * bs_coupling : Only Bending and stretching coupling
                            * bt_coupling : Only Bending and Twisting coupling
                            * bb_coupling : Only bending-1 and bending-2 coupling
                            * bend        : Both bending motions with their coupling
                            * st          : Stretching and twisting motions with their coupling
                            * bs          : Bending (b1, b2) and stretching motions with their coupling
                            * bt          : Bending (b1, b2) and twisting motions with their coupling

                            In case of elasticity type "ST", only following four energy terms are available "all", "diag", "stretch", "twist" and
                            "st_coupling".

                            The terms should provided as comma separated values. e.g. -et "full,diag,b1,b2,stretch,twist".

      -tbp total-bp-number, --total-bp total-bp-number
                            Total number of basepair in DNA/RNA.
                            It is an essential input.

      -estype esType, --elasticity-type esType
                            Elastic Properties type.
                            Two keywords are accepted: "BST" or "ST".
                            * "BST" : Bending-Stretching-Twisting --- All motions are considered
                            " "ST"  : Stretching-Twisting --- Bending motions are ignored.

                            WARNING: For accurate calculation of bending motions, DNA structures in trajectory must
                            be superimposed on a reference structure (See Publication's Method Section).

      -bs bp/s-start-number, --bp-start bp/s-start-number
                            First BP/BPS of DNA after which parameter will be extracted.
                            If it is not given, first basepair or base-step will be considered.

      -be bp/s-end-number, --bp-end bp/s-end-number
                            Last BP/BPS of DNA upto which parameter will be extracted.

                            If it is not given, last basepair or base-step will be considered.

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

See example `here <../global_elasticity.html#Global-deformation-free-energy>`_.