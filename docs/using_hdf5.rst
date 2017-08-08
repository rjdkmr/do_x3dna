Speed up dnaMD with HDF5 file
=============================

dnaMD tools and Python module need do_x3dna output files for the further analysis and calculation.
However, these do_x3dna output files could be large due to long trajectory. Reading
these large files could take **long time** and may become bottleneck during analysis. Moreover,
when these files are read, all data are kept in memory (RAM), and for long trajectory, dnaMD
may **consume considerable RAM**.

To overcome these limitations, we have implemented to use a file where all data can be stored for further
analysis and calculation. The file is in `HDF format <https://support.hdfgroup.org/HDF5/>`_.

Advantages of using HDF file
----------------------------

* Very fast to read and write.
* Data is indexed --- only required data can be read when neccessary.
* Can be read with other programming languages.
* Data are stored in disk --- reduce RAM consumption.

As can be seen below in examples, we no longer require the do_x3dna output files because 
all data are stored in HDF5 file. Therefore, data can be rapidly extracted from this HDF5 
file and can be further processed according to the requirement. Since, only few parameters
are neccessary at once during analysis, only these few parameters can be extracted
to reduce the memory consumption.


Using HDF5 file with dnaMD Python module
----------------------------------------
To store the data in HDF5 file, just initialize the class with the filename as shown below. 
Here, we named the HDF5 file as ``dna.h5``.

.. code-block:: Python

    dna = dnaMD.DNA(60, filename='dna.h5')     #Initialization for 60 base-pairs DNA bound with the protein


The usage of all methods and functions will be same because the file will be automatically used 
as a container to store and to retrieve the parameter values. 

The detail has been described in `this tutorial <notebooks/hdf5_tutorial.html>`_.

Using HDF5 file with dnaMD tools
--------------------------------

`saveAsH5 <commands/saveAsH5.html>`_ can be used to store the data from do_x3dna output files.
For example, to store the data from all do_x3dna ouput files to ``cdna.h5`` file, following command
can be used.

.. code-block:: bash

    dnaMD saveAsH5  -tbp 60 \ 
                    -lbp tutorial_data/L-BP_cdna.dat  \
                    -lbps tutorial_data/L-BPS_cdna.dat \
                    -lbph tutorial_data/L-BPH_cdna.dat \
                    -bbd tutorial_data/BackBoneCHiDihedrals_cdna.dat \
                    -mg tutorial_data/MGroove_cdna.dat \
                    -ha tutorial_data/HelAxis_cdna.dat \
                    -hr tutorial_data/HelixRad_cdna.dat \
                    -o cdna.h5


Subsequently, global helical axis, curvature and tangent can be calculated with `axisCurv <commands/axisCurv.html>`_ command.

.. code-block:: bash

    dnaMD axisCurv -tbp 60 -bs 4 -be 56 -io cdna.h5 -ctan -ap cdna_axis.pdb -scp 1.0

To plot the curvature values between 25th to 35th base-pairs, curvature as a function of time
can be extracted from HDF5 file by `vsTime <commands/vsTime.html>`_.

.. code-block:: bash

    dnaMD vsTime -i cdna.h5 -tbp 60 -bs 25 -be 35 -p "curvature" -mm sum -o curv.dat

Following output is obtained in ``curv.dat`` file.

::

    # Time 	 "curvature"
    0.0	    0.048152036382336755
    100.0	0.04615203638233675
    200.0	0.07015203638233675
    300.0	0.07415203638233675
    400.0	0.08415203638233676
    500.0	0.07515203638233674
    600.0	0.049152036382336756
    700.0	0.06115203638233675
    800.0	0.04613147157580158
    ...
    ...
    ...


