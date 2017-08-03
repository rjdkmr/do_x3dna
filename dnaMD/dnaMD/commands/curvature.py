import os
import sys
import re
import argparse

import dnaMD

description=\
"""Clculate global helical-axis, curvatures and tangents from local helical axis.
=================================================================================
It can be used to calulate global helical-axis and further curvatures and
tangents along it.

The local helical axis coordinates should be present in input HDF5 file. It can
be stored in HDF5 file using 'dnaMD saveasH5' command with '-ha/--helical-axis'
option.

At first global helical axis is calulated by spline interpolation using local
helical axis. This method smoothen the local helical axis and can be used as
global helical axis.

After determining the global helical axis, "-ctan/--curv-tangent" option can be
used to calculate curvature and tangents at each helical axis position
(each base-step). Curvature and tangents data will be stored in same HDF5 file.

This HDF5 file can be used later to extract and calculate data such as bending
angle, overall curvature etc.

"""

totalBpHelp=\
"""Total number of basepair in DNA/RNA.
It is an essential input.

"""

bpFirstHelp=\
"""Basepair number of first base-pair.
Usually it is one. Therefore, if this option is not provided, base-pair
numbering will start from one.

In rare cases, base-pair numbering might start with other number. In those
cases, use this option to start numbering of basepair from other number than
one.

"""

bpStartHelp=\
"""First BP/BPS of considered DNA segment.

If it is not given, first base-step will be considered.

"""

bpEndHelp=\
"""Last BP/BPS of considered DNA segment.

If it is not given, last base-step will be considered.

"""

inoutFileHelp=\
""" Name of input-output file.
It is a HDF5 format file generated from saveAsH5 tool. It should contain local
helical axis coordinates generated from do_x3dna.

It is also the output file as global helical axis, curvature and tangents will
be stored/saved in same file.

"""

curvTangentHelp=\
""" Calculate curvature and tangents.

If this options is used, curvature and tangents will be calculated from
global helical axis and stored in same HDF5 file.

If this option is not used, curvature and tangents will not be calculated.

The tangets can be used later to calculate bending angle.

"""

smoothHelp=\
""" A smoothing condition for spline interpolation.
For more details, see below in the link about "s = None", which is paased into
"scipy.interpolate.splprep()" function:
http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.splprep.html#scipy.interpolate.splprep

NOTE:
=> Lower value may lead to an artifact of local sharp kink in the smoothed axis.
=> Higher value may lead to the calculation of wrong helical axis.

"""

splineHelp=\
"""Degree of spline.
For more details, see below in the link about "k = 3", which is paased into
"scipy.interpolate.splprep()" function:
http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.splprep.html#scipy.interpolate.splprep

Its value should be from one to five.

"""

fillpointHelp=\
"""Number of interpolated points between two adjacent helical-axis coordinates.
Its value should be from one to ten.

"""

cutoffAngleHelp=\
"""Cut-off bending angle to define sharp kink in fitted 3D curve or axis.
If angle in fitted curve is larger than this cut-off, refitting will be
performed after deleting few of the original helical axis positions.
If after this deletions, bending angle will not reduce below cut-off angle,
value of ""smooth"" will be increased by 100 and entire cycle of
fitting-refitting will be performed. When, value of ""smooth"" increases to
more than 10000 during this fitting-refitting cycles, fitting process will be
stopped with a warning message.

If fitting process is not successfull (angle > cutoff-angle), a mask will be
added for the respective frame. This mask can be used later to discard these
frames during subsequent analysis.

"""

axisPdbHelp=\
""" Output helical-axis as PDB file.

This option enable the output of helical axis to a PDB file. Either local or
global or both axis can be included in output file. Also, scaled curvature can
be included in b-factor field of PDB file.

"""

axisPdbOptionHelp=\
""" Option for helical-axis PDB file.

Following three options can be used to modify the output in PDB file:
* 'global': Only global smoothed helical axis will be presen.
* 'local' : Only local helical axis calculated from do_x3dna will be present.
* 'both' : Both global and local helical axis will be present.

"""

def main():
    parser, args = parseArguments()

    # Total number of base-pair
    firstBP = args.firstBP
    totalBP = None
    if args.totalBP is None:
        showErrorAndExit(parser, "No total number of BP!!!\n")
    else:
        totalBP = args.totalBP

    toMinusBP=2   # Base-step
    # Determine start and end-bp
    startBP = None
    if args.startBP is None:
        startBP = args.firstBP
    else:
        startBP = args.startBP

    endBP = None
    if  args.endBP is None:
        endBP = totalBP+firstBP-toMinusBP
    else:
        endBP = args.endBP

    # Determine step_range and step
    step_range = False
    step = None
    if not (startBP == args.firstBP and endBP == totalBP+firstBP-toMinusBP):
        # Check consistency of start bp
        if (startBP < firstBP) or (startBP > totalBP+firstBP-toMinusBP):
            msg = 'The requested start bp {0} is out side of {1}-{2} range.'.format(startBP, firstBP, totalBP+firstBP-toMinusBP)
            showErrorAndExit(parser, msg)

        # Check consistency of end-bp
        if startBP > endBP:
            msg = 'The requested end bp {0} is larger than requested start bp {1}!!!'.format(endBP, startBP)
            showErrorAndExit(parser, msg)

        if (endBP > totalBP+firstBP-toMinusBP) or (endBP < firstBP):
            msg = 'The requested end bp {0} is out side of {1}-{2} range.'.format(endBP, firstBP, totalBP+firstBP-toMinusBP)
            showErrorAndExit(parser, msg)

        step_range = True
        step = [startBP, endBP]

    # Determine file-extension type
    fileType = 'hdf5'
    ouputFileExtension = os.path.splitext(args.ioFile)
    if ouputFileExtension not in ['.h5', '.hdf5', 'hdf']:
        showErrorAndExit(parser, "File extension {0} is not recognized as an \
        aceeptable HDF5 extension.\n Use '.h5', '.hdf5' \
        or '.hdf'.".format(ouputFileExtension))

    # initialize DNA object
    dna = dnaMD.DNA(totalBP, filename=args.ioFile, startBP=firstBP)
    dna.generate_smooth_axis(step_range=step_range, step=step,
                            smooth=args.smooth,
                            spline=args.spline, fill_point=args.fill_point,
                            cut_off_angle=args.cut_off_angle)

    if args.curvTangent:
        print(" ... Calculating curvature and tangents ")
        dna.calculate_curvature_tangent(step_range=step_range, step=step,
                                        store_tangent=True)
        print("                              ... Finished")

def parseArguments():
    parser = argparse.ArgumentParser(prog='dnaMD smoothenAxis',
                                    description=description,
                                    formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-tbp', '--total-bp', action='store',
                        type=int, metavar='total-bp-number',
                        dest='totalBP', help=totalBpHelp)

    parser.add_argument('-fbp', '--first-bp', action='store',
                        type=int, metavar='1', default=1,
                        dest='firstBP', help=bpFirstHelp)

    parser.add_argument('-io', '--in-out-put', action='store',
                        type=str, metavar='inout.h5',
                        dest='ioFile', help=inoutFileHelp)

    parser.add_argument('-bs', '--bp-start', action='store',
                        type=int, metavar='bp/s-start-number',
                        default=1,
                        dest='startBP', help=bpStartHelp)

    parser.add_argument('-be', '--bp-end', action='store',
                        type=int, dest='endBP',
                        metavar='b/ps-end-number',
                        help=bpEndHelp)

    parser.add_argument('-ctan', '--curv-tangent', action='store_true',
                        dest='curvTangent', help=curvTangentHelp)

    parser.add_argument('-s', '--smooth', action='store',
                        type=float, dest='smooth',
                        metavar='500.0',default=500.0,
                        help=smoothHelp)

    parser.add_argument('-st', '--smooth', action='store',
                        type=float, dest='smooth',
                        metavar='500.0',default=500.0,
                        help=smoothHelp)

    parser.add_argument('-sp', '--spline', action='store',
                        type=int, dest='spline',
                        metavar=3,default=3,choices=[1,2,3,4,5],
                        help=splineHelp)

    parser.add_argument('-fp', '--fill-point', action='store',
                        type=int, dest='fill_point',
                        metavar=6,default=6, choices=range(1,stop=10),
                        help=fillpointHelp)

    parser.add_argument('-cta', '--cut-off-angle', action='store',
                        type=float, dest='cut_off_angle',
                        metavar=20,default=20,
                        help=cutoffAngleHelp)

    parser.add_argument('-ap', '--axis-pdb', action='store',
                        type=float, dest='pdbFile',
                        metavar='helical_axis.pdb',default='helical_axis.pdb',
                        help=axisPdbHelp)

    parser.add_argument('-apo', '--axis-pdb-option', action='store',
                        type=float, dest='pdbAxisOption',
                        metavar='both',default='both',
                        choices=['global', 'local', 'both'],
                        help=axisPdbOptionHelp)


    idx = sys.argv.index("smoothenAxis")+1
    args = parser.parse_args(args=sys.argv[idx:])

    return parser, args


def showErrorAndExit(parser, message):
    parser.print_help()
    print("\n===== ERROR =======")
    print(message)
    print("See Usage Above!!!")
    sys.exit(False)
