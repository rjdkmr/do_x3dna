#!/usr/bin/env python
#
#
# This file is part of do_x3dna
#
# Author: Rajendra Kumar
# Copyright (C) 2014-2017  Rajendra Kumar
#
# do_x3dna uses 3DNA package (http://x3dna.org).
# Please cite the original publication of the 3DNA package:
# Xiang-Jun Lu & Wilma K. Olson (2003)
# 3DNA: a software package for the analysis, rebuilding and visualization of
# three-dimensional nucleic acid structures
# Nucleic Acids Res. 31(17), 5108-21.
#
# do_x3dna is a free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# do_x3dna is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with do_x3dna.  If not, see <http://www.gnu.org/licenses/>.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
# TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#============================================================================

import os
import sys
import re
import argparse

import dnaMD

description=\
"""Calculate global helical-axis, curvatures and tangents from local helical axis.
=================================================================================
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

The tangents can be used later to calculate bending angle.

"""

smoothHelp=\
""" A smoothing condition for spline interpolation.
For more details, see below in the link about "s = None", which is passed into
"scipy.interpolate.splprep()" function:
http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.splprep.html#scipy.interpolate.splprep

NOTE:
=> Lower value may lead to an artifact of local sharp kink in the smoothed axis.
=> Higher value may lead to the calculation of wrong helical axis.

"""

splineHelp=\
"""Degree of spline.
For more details, see below in the link about "k = 3", which is passed into
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

If fitting process is not successful (angle > cutoff-angle), a mask will be
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
* 'global': Only global smoothed helical axis will be present.
* 'local' : Only local helical axis calculated from do_x3dna will be present.
* 'both' : Both global and local helical axis will be present.

"""

scaleCurvPdbHelp=\
""" Store curvature in PDB file

This option writes the curvature values of smoothed global helical
axis in B-factor column of PDB file. This value further scale the
calculated curvature value in PDB file.

For example: '-scp 10' option will multiply all curvature values with
10 and write it to PDB file's B-factor column.

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
    outputFileExtension = os.path.splitext(args.ioFile)[1]
    if outputFileExtension not in ['.h5', '.hdf5', 'hdf']:
        showErrorAndExit(parser, "File extension {0} is not recognized as an \
        acceptable HDF5 extension.\n Use '.h5', '.hdf5' \
        or '.hdf'.".format(outputFileExtension))

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

    if args.pdbFile is not None:

        if args.pdbAxisOption == 'both':
            write_smooth_axis = True
            write_orig_axis = True
        elif pdbAxisOption == 'global':
            write_smooth_axis = True
            write_orig_axis = False
        elif pdbAxisOption == 'local':
            write_smooth_axis = False
            write_orig_axis = True
        else:
            write_smooth_axis = True
            write_orig_axis = False

        if args.scaleCurvPdb is None:
            write_curv = False
            args.scaleCurvPdb = 1
        else:
            write_curv = True

        dna.write_haxis_pdb(filename=args.pdbFile,
                            step_range=step_range,
                            step=step, write_smooth_axis=write_smooth_axis,
                            write_orig_axis=write_orig_axis,
                            write_curv=write_curv,
                            scale_curv=args.scaleCurvPdb)

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

    parser.add_argument('-sp', '--spline', action='store',
                        type=int, dest='spline',
                        metavar=3,default=3,choices=[1,2,3,4,5],
                        help=splineHelp)

    parser.add_argument('-fp', '--fill-point', action='store',
                        type=int, dest='fill_point',
                        metavar=6,default=6, choices=list(range(1,10)),
                        help=fillpointHelp)

    parser.add_argument('-cta', '--cut-off-angle', action='store',
                        type=float, dest='cut_off_angle',
                        metavar=20,default=20,
                        help=cutoffAngleHelp)

    parser.add_argument('-ap', '--axis-pdb', action='store',
                        type=str, dest='pdbFile',
                        metavar='helical_axis.pdb',default='helical_axis.pdb',
                        help=axisPdbHelp)

    parser.add_argument('-apo', '--axis-pdb-option', action='store',
                        type=str, dest='pdbAxisOption',
                        metavar='both',default='both',
                        choices=['global', 'local', 'both'],
                        help=axisPdbOptionHelp)

    parser.add_argument('-scp', '--scale-curv-pdb', action='store',
                        type=float, dest='scaleCurvPdb',
                        metavar=1.0,default=None,
                        help=scaleCurvPdbHelp)


    idx = sys.argv.index("axisCurv")+1
    args = parser.parse_args(args=sys.argv[idx:])

    return parser, args


def showErrorAndExit(parser, message):
    parser.print_help()
    print("\n===== ERROR =======")
    print(message)
    print("See Usage Above!!!")
    sys.exit(False)
