axisCurv
========

**Calculate global helical-axis, curvatures and tangents from local helical axis.**

It can be used to calculate global helical-axis and further curvatures and
tangents along it.

The local helical axis coordinates should be present in input HDF5 file. It can
be stored in HDF5 file using 'dnaMD saveasH5' command with '-ha/--helical-axis'
option.

At first global helical axis is calculated by spline interpolation using local
helical axis. This method makes the local helical axis smooth and can be used as
global helical axis.

After determining the global helical axis, "-ctan/--curv-tangent" option can be
used to calculate curvature and tangents at each helical axis position
(each base-step). Curvature and tangents data will be stored in same HDF5 file.

This HDF5 file can be used later to extract and calculate data such as bending
angle, overall curvature etc.


**Usage:**

::

  usage: dnaMD smoothenAxis [-h] [-tbp total-bp-number] [-fbp 1] [-io inout.h5]
                            [-bs bp/s-start-number] [-be b/ps-end-number]
                            [-ctan] [-s 500.0] [-sp 3] [-fp 6] [-cta 20]
                            [-ap helical_axis.pdb] [-apo both] [-scp 1.0]


**Optional arguments:**

::

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
                          
    -io inout.h5, --in-out-put inout.h5
                          Name of input-output file.
                          It is a HDF5 format file generated from saveAsH5 tool. It should contain local
                          helical axis coordinates generated from do_x3dna.
                          
                          It is also the output file as global helical axis, curvature and tangents will
                          be stored/saved in same file.
                          
    -bs bp/s-start-number, --bp-start bp/s-start-number
                          First BP/BPS of considered DNA segment.
                          
                          If it is not given, first base-step will be considered.
                          
    -be b/ps-end-number, --bp-end b/ps-end-number
                          Last BP/BPS of considered DNA segment.
                          
                          If it is not given, last base-step will be considered.
                          
    -ctan, --curv-tangent
                          Calculate curvature and tangents.
                          
                          If this options is used, curvature and tangents will be calculated from
                          global helical axis and stored in same HDF5 file.
                          
                          If this option is not used, curvature and tangents will not be calculated.
                          
                          The tangents can be used later to calculate bending angle.
                          
    -s 500.0, --smooth 500.0
                          A smoothing condition for spline interpolation.
                          For more details, see below in the link about "s = None", which is passed into
                          "scipy.interpolate.splprep()" function:
                          http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.splprep.html#scipy.interpolate.splprep
                          
                          NOTE:
                          => Lower value may lead to an artifact of local sharp kink in the smoothed axis.
                          => Higher value may lead to the calculation of wrong helical axis.
                          
    -sp 3, --spline 3     Degree of spline.
                          For more details, see below in the link about "k = 3", which is passed into
                          "scipy.interpolate.splprep()" function:
                          http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.splprep.html#scipy.interpolate.splprep
                          
                          Its value should be from one to five.
                          
    -fp 6, --fill-point 6
                          Number of interpolated points between two adjacent helical-axis coordinates.
                          Its value should be from one to ten.
                          
    -cta 20, --cut-off-angle 20
                          Cut-off bending angle to define sharp kink in fitted 3D curve or axis.
                          If angle in fitted curve is larger than this cut-off, refitting will be
                          performed after deleting few of the original helical axis positions.
                          If after this deletions, bending angle will not reduce below cut-off angle,
                          value of ""smooth"" will be increased by 100 and entire cycle of
                          fitting-refitting will be performed. When, value of ""smooth"" increases to
                          more than 10000 during this fitting-refitting cycles, fitting process will be
                          stopped with a warning message.
                          
                          If fitting process is not successful (angle > cutoff-angle), a mask will be
                          added for the respective frame. This mask can be used later to discard these
                          frames during subsequent analysis.
                          
    -ap helical_axis.pdb, --axis-pdb helical_axis.pdb
                          Output helical-axis as PDB file.
                          
                          This option enable the output of helical axis to a PDB file. Either local or
                          global or both axis can be included in output file. Also, scaled curvature can
                          be included in b-factor field of PDB file.
                          
    -apo both, --axis-pdb-option both
                          Option for helical-axis PDB file.
                          
                          Following three options can be used to modify the output in PDB file:
                          * 'global': Only global smoothed helical axis will be present.
                          * 'local' : Only local helical axis calculated from do_x3dna will be present.
                          * 'both' : Both global and local helical axis will be present.
                          
    -scp 1.0, --scale-curv-pdb 1.0
                          Store curvature in PDB file
                          
                          This option writes the curvature values of smoothed global helical 
                          axis in B-factor column of PDB file. This value further scale the
                          calculated curvature value in PDB file.
                          
                          For example: '-scp 10' option will multiply all curvature values with
                          10 and write it to PDB file's B-factor column.
  