#!/usr/bin/env python
#
#
# This file is part of do_x3dna
#
# Author: Rajendra Kumar
# Copyright (C) 2014-2023  Rajendra Kumar
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

description="""DESCRIPTION
===========
Global Elastic Properties of the DNA

This can be used to Global Elastic Properties of the DNA from the simulations.

"""

inputFileHelp=\
"""Name of input file (hdf5 file).
This file should contain the required parameters. It should be hdf5 storage file.

"""

outputFileHelp=\
"""Name of output file with csv extension.
This file will contain the elasticity modulus matrix where values will be separated 
by comma. Since modulus matrix is also shown as screen output, this option is not 
necessary.

"""

vsTimeHelp=\
"""Calculate elasticity as a function of time and save in this csv format file.
It can be used to obtained elastic moduli as a function of time to check their 
convergence. If this option is used, -fgap/--frame-gap is an essential option.

NOTE: Elastic properties cannot be calculated using a single frame because 
fluctuations are required. Therefore, here time means trajectory between zero 
time to given time.

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
esTypeHelp=\
"""Elastic Properties type.
Two keywords are accepted: "BST" or "ST".
* "BST" : Bending-Stretching-Twisting --- All motions are considered
* "ST"  : Stretching-Twisting --- Bending motions are ignored.

WARNING: For accurate calculation of bending motions, DNA structures in trajectory must 
be superimposed on a reference structure (See Publication's Method Section).

"""

bpStartHelp=\
"""First BP/BPS of DNA after which parameter will be extracted.
If it is not given, first basepair or base-step will be considered.

"""

bpEndHelp=\
"""Last BP/BPS of DNA upto which parameter will be extracted.

If it is not given, last basepair or base-step will be considered.

"""

frameGapHelp=\
    """Number of frames to skip for next calculation
When calculating elastic modulus as a function of time, this option will determine
the time-gap between each calculation point. Lower the time-gap, slower will be the 
calculation.
"""

paxisHelp=\
"""Principal axis parallel to global helical-axis
Three keywords are accepted: "X", "Y" and "Z". Only require when bending motions are
included in the calculation.

"""

errorHelp=\
""" Error of elastic modulus
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

"""

toolsHelp=\
"""Tools to calculate autocorrelation time or bloack averaging error.
By default it is g_analyze (Gromacs-4.5.x/4.6.x versions). For newer versions, use "gmx analyze".

"""

def main():
    parser, args = parseArguments()

    # Input file
    inputFile = None
    if args.inputFile is not None:
        inputFile = args.inputFile.name
        args.inputFile.close()
    else:
        showErrorAndExit(parser, "No Input File!!!\n")

    # Determine file-extension type
    inputFileExtension = os.path.splitext(inputFile)[1]
    if inputFileExtension not in ['.h5', '.hdf5', '.hdf']:
        showErrorAndExit(parser, "File should be in HDF5 (h5, hdf5 or hdf file extension) format...\n")


    # Total number of base-pair
    firstBP = args.firstBP
    totalBP = None
    if args.totalBP is None:
        showErrorAndExit(parser, "No total number of BP!!!\n")
    else:
        totalBP = args.totalBP

    # Determine start and end-bp
    toMinusBP = 2
    startBP = args.startBP
    if startBP is None:
        startBP = firstBP
    endBP = args.endBP
    if endBP is None:
        endBP = firstBP + totalBP - toMinusBP

    # Check consistency of start bp
    if (startBP < firstBP) or (startBP > totalBP+firstBP-toMinusBP):
        msg = 'The requested start bp {0} is out side of {1}-{2} range.'.format(startBP, firstBP, totalBP+firstBP-toMinusBP)
        showErrorAndExit(parser, msg)

    # Check consistency of end-bp
    if endBP is not None:
        if startBP > endBP:
            msg = 'The requested end bp {0} is larger than requested start bp {1}!!!'.format(endBP, startBP)
            showErrorAndExit(parser, msg)

        if (endBP > totalBP+firstBP-toMinusBP) or (endBP < firstBP):
            msg = 'The requested end bp {0} is out side of {1}-{2} range.'.format(endBP, firstBP, totalBP+firstBP-toMinusBP)
            showErrorAndExit(parser, msg)

    # Define DNA segement here
    bp = [startBP, endBP]

    if args.esType == 'BST' and args.paxis is None:
        showErrorAndExit(parser, 'To calculate bending, principal axis parallel to helical axis is required.')

    if args.vsTime is not None and args.frameGap is None:
        showErrorAndExit(parser, 'Frame-gap is required when elasticity as a function of time is calculated.')

    if args.err_type is not None and args.frameGap is None:
        showErrorAndExit(parser, 'Frame-gap is required when error is calculated.')

    # initialize DNA object
    dna = dnaMD.dnaEY(totalBP, esType=args.esType, filename=inputFile, startBP=firstBP)

    # Check if mask is in object
    if dna.dna.mask is not None:
        masked = True
    else:
        masked = False

    if args.vsTime is None:
        out = ''
        towrite = ''
        mean_out = output_props(args.esType)
        if args.esType == 'BST':
            mean, result = dna.getStretchTwistBendModulus(bp, paxis=args.paxis, masked=True, matrix=False)
            for i in range(4):
                for j in range(4):
                    if j != 3:
                        out += '{0:>10.3f}  '.format(result[i][j])
                        towrite += '{0:.3f},'.format(result[i][j])
                    else:
                        out += '{0:>10.3f}\n'.format(result[i][j])
                        towrite += '{0:.3f}\n'.format(result[i][j])

                mean_out += '{0:>15.3f}  '.format(mean[i])
            mean_out += '\n'

        if args.esType == 'ST':
            mean, result = dna.getStretchTwistModulus(bp, masked=masked, matrix=False)
            for i in range(2):
                for j in range(2):
                    if j != 1:
                        out += '{0:.3f}  '.format(result[i][j])
                        towrite += '{0:.3f},'.format(result[i][j])
                    else:
                        out += '{0:.3f}\n'.format(result[i][j])
                        towrite += '{0:.3f}\n'.format(result[i][j])

                mean_out += '{0:>12.3f}  '.format(mean[i])
            mean_out += '\n'

        print('=========== Elastic Modulus Matrix ===========\n')
        print(out)
        print('=========== ====================== ===========\n')

        print('=========================  Average Values  ==========================\n')
        print(mean_out)
        print('=========== ====================== ====================== ===========')

        if args.outputFile is not None:
            with open(args.outputFile, 'w') as fout:
                fout.write(towrite)

    else:
        time, modulus = dna.getModulusByTime(bp, args.frameGap, masked=masked, paxis=args.paxis, outFile=args.vsTime)

        props_name = list(modulus.keys())
        if args.err_type is not None:
            error = dnaMD.get_error(time, list(modulus.values()), len(props_name), err_type=args.err_type, tool=args.tool)
            sys.stdout.write("==============================================\n")
            sys.stdout.write('{0:<16}{1:>14}{2:>14}\n'.format('Elasticity', 'Value', 'Error'))
            sys.stdout.write("----------------------------------------------\n")
            for i in range(len(props_name)):
                sys.stdout.write('{0:<16}{1:>14.3f}{2:>14.3f}\n'.format(props_name[i], modulus[props_name[i]][-1], error[i]))
            sys.stdout.write("==============================================\n\n")


def output_props(key):
    out = ''
    if key == 'BST':
        props = ['Bending-1 Angle', 'Bending-2 Angle', 'Contour Length', 'Sum. Twist']
        for i in range(4):
            out += props[i] + '    '
        out += '\n'

    if key == 'ST':
        props = ['Contour Length', 'Sum. Twist']
        for i in range(2):
            out += props[i] + '   '
        out += '\n'

    return out


def parseArguments():
    parser = argparse.ArgumentParser(prog='dnaMD globalElasticity',
                                     description=description,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', action='store',
                        type=argparse.FileType('rb'), metavar='parameter.h5',
                        dest='inputFile', required=False, help=inputFileHelp)

    parser.add_argument('-o', '--output', action='store',
                        type=str, metavar='output.csv',
                        dest='outputFile', help=outputFileHelp)

    parser.add_argument('-ot', '--output-time', action='store',
                        type=str, metavar='elasicity_vs_time.csv',
                        dest='vsTime', help=vsTimeHelp)

    parser.add_argument('-tbp', '--total-bp', action='store',
                        type=int, metavar='total-bp-number',
                        dest='totalBP', help=totalBpHelp)

    parser.add_argument('-estype', '--elasticity-type', action='store',
                        dest='esType', metavar='esType', default='ST',
                        choices = ['ST', 'BST'],
                        type=str, help=esTypeHelp)

    parser.add_argument('-bs', '--bp-start', action='store',
                        type=int, metavar='bp/s-start-number',
                        default=1,
                        dest='startBP', help=bpStartHelp)

    parser.add_argument('-be', '--bp-end', action='store',
                        type=int, dest='endBP',
                        metavar='bp/s-end-number',
                        help=bpEndHelp)

    parser.add_argument('-fgap', '--frame-gap', action='store',
                        type=int, dest='frameGap',
                        metavar='frames-to-skip',
                        help=frameGapHelp)

    parser.add_argument('-paxis', '--principal-axis', action='store',
                        type=str, metavar='X', default=None,
                        choices=['X', 'Y', 'Z'],
                        dest='paxis', help=paxisHelp)

    parser.add_argument('-em', '--error-method', action='store',
                        type=str, metavar='block', default=None,
                        choices=['std', 'acf', 'block'],
                        dest='err_type', help=errorHelp)

    parser.add_argument('-gt', '--gromacs-tool', action='store',
                        type=str, metavar='gmx analyze', default='gmx analyze',
                        dest='tool', help=toolsHelp)

    parser.add_argument('-fbp', '--first-bp', action='store',
                        type=int, metavar='1', default=1,
                        dest='firstBP', help=bpFirstHelp)

    idx = sys.argv.index("globalElasticity") + 1
    args = parser.parse_args(args=sys.argv[idx:])

    return parser, args

def showErrorAndExit(parser, message):
    parser.print_help()
    print("\n===== ERROR =======")
    print(message)
    print("See Usage Above!!!")
    sys.exit(False)

if __name__ == '__main__':
    main()
