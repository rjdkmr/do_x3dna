saveAsH5
========

**Save parameters to a HDF5 file**

The parameters can be read from the do_x3dna output files and stored in a HDF5
file. This step enables rapid reading and processing for subsequent analysis of
data.

Several calculations such as determining global helical axis, its curvature and
bending angles needs this HDF5 file because helical axis, curvature and
tangents data will be stored in this HDF5 file.

**Usage:**

.. code-block:: bash

    usage: dnaMD saveAsH5 [-h] [-tbp total-bp-number] [-fbp 1]
                          [-i input1.dat,input2.dat,input3.dat] [-o output.h5]


**optional arguments:**

.. code-block:: bash

      -h, --help            show this help message and exit
      -tbp total-bp-number, --total-bp total-bp-number
                            Total number of basepair in DNA/RNA.
                            It is an essential input.

      -fbp 1, --first-bp 1  Basepair number of first base-pair.
                            Usually it is one. Therefore, if this option is not provided, base-pair
                            numbering will start from one.

                            In rare cases, base-pair numbering might start with other number. In those
                            cases, use this option to start numbering of basepair from other number than
                            one.

      -i input1.dat,input2.dat,input3.dat, --input-files input1.dat,input2.dat,input3.dat
                            List of input files containing parameters.
                            List of input files separated by comma (,) containing parameters data.
                            These file should be a output file from do_x3dna.
                            Following files are accepted:
                            => L-BP_*.dat
                            => L-BPS_*.dat
                            => L-BPH_*.dat
                            => HelAxis_*.dat
                            => MGroove_*.dat
                            => HelixRad_*.dat
                            => BackBoneCHiDihedrals_*.dat

      -o output.h5, --output output.h5
                             Name of output file.
                            It is a hdf5 format file. Therefore, use ".h5" extension with it.


Example
-------

`saveAsH5` accepts list of input files separated by comma (,) as following:

.. code-block:: bash

    dnaMD saveAsH5  -tbp 60 \
                    -i tutorial_data/L-BP_cdna.dat,tutorial_data/L-BPS_cdna.dat,tutorial_data/L-BPH_cdna.dat,\
                    tutorial_data/BackBoneCHiDihedrals_cdna.dat,tutorial_data/MGroove_cdna.dat,\
                    tutorial_data/HelAxis_cdna.dat,tutorial_data/HelixRad_cdna.dat \
                    -o pdna.h5



Above command load all the parameters from these files and store in the output HDF5 file.