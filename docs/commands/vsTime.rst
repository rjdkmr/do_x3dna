vsTime
======

**Parameter as a function of time**

This can be used to extract the parameter of either a specific base-pair/step
or over a DNA segment as a function of time.

**Usage:**

.. code-block:: bash

    usage: dnaMD vsTime [-h] [-i L-BP_cdna.dat] [-o output.dat]
                        [-tbp total-bp-number] [-fbp 1] [-p parameter]
                        [-bs bp/s-start-number] [-be b/ps-end-number]
                        [-mm sum-or-mean]


**Optional arguments:**

.. code-block:: bash

    -h, --help            show this help message and exit
    -i L-BP_cdna.dat, --input L-BP_cdna.dat
                          Name of input file.
                          This file should contain the required parameters.

    -o output.dat, --output output.dat
                          Name of ouput file.
                          The extracted output will be written in output file.

    -tbp total-bp-number, --total-bp total-bp-number
                          Total number of basepair in DNA/RNA.
                          It is an essential input.

    -fbp 1, --first-bp 1  Basepair number of first base-pair.
                          Usually it is one. Therefore, if this option is not provided, base-pair
                          numbering will start from one.

                          In rare cases, base-pair numbering might start with other number. In those
                          cases, use this option to start numbering of basepair from other number than
                          one.

    -p parameter, --parameter parameter
                          Parameter name.
                          This paramter will be extracted from file. Ensure that parameter is present
                          in the file, otherwise wrong values will be extracted from file.

    -bs bp/s-start-number, --bp-start bp/s-start-number
                          First BP/BPS of a segment for extraction
                          If it is not given, first basepair or base-step will be considered.

    -be b/ps-end-number, --bp-end b/ps-end-number
                          Last BP/BPS of a segment for extraction

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
