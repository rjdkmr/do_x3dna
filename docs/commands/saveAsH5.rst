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

    usage: dnaMD saveAsH5 [-h] [-tbp total-bp-number] [-fbp 1] [-lbp L-BP_g.dat]
                          [-lbps L-BPS_g.dat] [-lbph L-BPH_g.dat]
                          [-bbd BackBoneCHiDihedrals_g.dat] [-mg MGroove_g.dat]
                          [-ha HelAxis_g.dat] [-hr HelixRad_g.dat] [-o output.h5]


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

    -lbp L-BP_g.dat, --local-base-pair L-BP_g.dat
                          File containing base-pair parameters.
                          This file is obtained from do_x3dna. It should contain six parameters:
                          (1) shear, (2) stretch, (3) stagger, (4) buckle, (5) propeller and (6) opening.

    -lbps L-BPS_g.dat, --local-base-step L-BPS_g.dat
                           File containing base-step parameters.
                          This file is obtained from do_x3dna. It should contain six parameters:
                          (1) shift, (2) slide, (3) rise, (4) tilt, (5) roll and (6) twist

    -lbph L-BPH_g.dat, --local-base-helical-step L-BPH_g.dat
                           File containing helical base-step parameters.
                          This file is obtained from do_x3dna. It should contain six parameters:
                          (1) x-disp, (2) y-disp, (3) h-rise, (4) inclination, (5) tip and (6) h-twist

    -bbd BackBoneCHiDihedrals_g.dat, --backbone-dihedrals BackBoneCHiDihedrals_g.dat
                           File containing backbone and chi dihedral angles.
                          This file is obtained from do_x3dna. It should contain 14 parameters:
                          (1) alpha S1, (2) beta S1, (3) gamma S1, (4) delta S1, (5) epsilon S1,
                          (6) zeta S1,  (7) chi S1,  (8) alpha S2, (9) beta S2,  (10) gamma S2,
                          (11) delta S2, (12) epsilon S2, (13)zeta S2, and (14) chi S2

                          where S1 and S2 is first and second strand respectively.

    -mg MGroove_g.dat, --major-minor-grooves MGroove_g.dat
                           File containing major and minor grooves width.
                          This file is obtained from do_x3dna. It should contain four parameters:
                          (1) minor groove, (2) minor groove refined,
                          (3) major groove, (4) major groove refined.

    -ha HelAxis_g.dat, --helical-axis HelAxis_g.dat
                           File containing local helical axis coordinates.
                          This file is obtained from do_x3dna. It should contain coordinates of local
                          helical axis.

                          Note that global and smooth helical axis can be calculated using same HDF5
                          file as input.

    -hr HelixRad_g.dat, --helical-radius HelixRad_g.dat
                          File containing helical radius.
                          This file is obtained from do_x3dna. It should contain local helical radius.

    -o output.h5, --output output.h5
                           Name of output file.
                          It is a hdf5 format file. Therefore, use ".h5" extension with it.
