#!/usr/bin/env python
#
#
# This file is part of do_x3dna
#
# Author: Rajendra Kumar
# Copyright (C) 2014-2018  Rajendra Kumar
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
Local Elastic Properties of the DNA

This can be used to calculate local elastic properties of the DNA from the simulations. Here 
local DNA segment referred to less than 5 base-pair (4 basepair-step)long.

WARNING: The option "-o/--output" and "-ot/--output-time" cannot be used with "-os/--output-segments" 
         simultaneously due to incompatible usage of "-bs/--bps-start" and "-be/--bps-end". Both 
         "-o/--output" and "-ot/--output-time" options can be used simultaneously.
         
WARNING: In case of "-o/--output" and "-ot/--output-time" difference between "-bs/--bps-start" 
         and "-be/--bps-end" should be less than 4 while "--span" should be less than 4 in 
         case of "-os/--output-segments".

"""

inputFileHelp=\
"""Name of input file (hdf5 file).
This file should contain the required parameters. It should be hdf5 storage file.

"""

outputFileHelp=\
"""Name of output file with csv extension.
This file will contain the local elasticity matrix where values will be separated by comma. 
This is also shown as screen output.

"""

outputFileSegmentsHelp=\
"""Calculate local elasticity of consecutive overlapped DNA segments and save in this csv format file.
It enables the calculation of local elastic properties of small overlapped DNA segments with 
error. This file will contain the local elastic properties of small overlapped DNA segments of 
length given by "-s/--span". As error will be calculated for each segment, '-fgap/--frame-gap' 
is required.

"""

vsTimeHelp=\
"""Calculate elasticity as a function of time and save in this csv format file.
It can be used to obtained local elasticity as a function of time to check their 
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

helicalHelp=\
"""Enable helical base-step parameters
By default, elasticity for base step-parameters (shift, slide, rise, 
tilt, roll, twist) are calculated. If this option is used, helical
base-step parameters (x-disp, y-disp, h-rise, inclination, tip, h-twist)
will be used for elasticity.

"""

bpsStartHelp=\
"""First basepair-step (bp-s) of DNA after which parameter will be extracted.
If it is not given, first bp-s will be considered.

In case of "-o/--output" and "-ot/--output-time" this is firs bp-s of
the segment while for "-os/--output-segments", it is first bp-s of
first overlapped segments.

"""

bpsEndHelp=\
"""Last basepair-step (bp-s) of DNA upto which parameter will be extracted.

In case of "-o/--output" and "-ot/--output-time" this is last bp-s of
the segment while for "-os/--output-segments", it is last bp-s of
last overlapped segments.

"""

spanHelp=\
"""Length of overlapping (local) DNA segments. 
It is essential when "-os/--output-segments" is used. It should not be larger than four.

"""

frameGapHelp=\
    """Number of frames to skip for next calculation
When calculating elastic modulus as a function of time, this option will determine
the time-gap between each calculation point. Lower the time-gap, slower will be the 
calculation.
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

    if args.outputFileSegments is not None:
        if args.vsTime is not None:
            showErrorAndExit(parser, "Both -ot/--output-time and -os/--output-segments cannot used simultaneously. \n")
        if args.outputFile is not None:
            showErrorAndExit(parser, "Both -o/--output and -os/--output-segments cannot used simultaneously. \n")


    # Total number of base-pair
    firstBP = args.firstBP
    totalBP = None
    if args.totalBP is None:
        showErrorAndExit(parser, "No total number of BP!!!\n")
    else:
        totalBP = args.totalBP

    # Determine start and end-bp
    toMinusBP = 2
    startBPS = args.startBPS
    if startBPS is None:
        showErrorAndExit(parser, "-bs/--bps-start is missing... \n")
    endBPS = args.endBPS
    if endBPS is None:
        showErrorAndExit(parser, "-be/--bps-end is missing... \n")

    # Check consistency of start bp
    if (startBPS < firstBP) or (startBPS > totalBP+firstBP-toMinusBP):
        msg = 'The requested start bp-s {0} is out side of {1}-{2} range.'.format(startBPS, firstBP, totalBP+firstBP-toMinusBP)
        showErrorAndExit(parser, msg)

    # Check consistency of end-bp
    if startBPS > endBPS:
        msg = 'The requested end bp-s {0} is larger than requested start bp-s {1}!!!'.format(endBPS, startBPS)
        showErrorAndExit(parser, msg)

    if (endBPS > totalBP+firstBP-toMinusBP) or (endBPS < firstBP):
        msg = 'The requested end bp-s {0} is out side of {1}-{2} range.'.format(endBPS, firstBP, totalBP+firstBP-toMinusBP)
        showErrorAndExit(parser, msg)

    if (args.vsTime is not None) or (args.outputFile is not None):
        if endBPS - startBPS > 4:
            showErrorAndExit(parser,
                             "Difference between '-bs/--bps-start' and '-be/--bps-end' should be less than 4. \n")

    # Define DNA segement here
    bp = [startBPS, endBPS]

    if args.vsTime is not None and args.frameGap is None:
        showErrorAndExit(parser, 'Frame-gap is required when elasticity as a function of time is calculated.\n')
    if args.outputFileSegments is not None and args.frameGap is None:
        showErrorAndExit(parser, 'Frame-gap is required when elasticities of consecutive DNA segments are calculated.\n')
    if args.outputFileSegments is not None and args.err_type is None:
        showErrorAndExit(parser, '-em/--error-method cannot be none when elasticities of consecutive DNA '
                                 'segments are calculated.\n')

    # initialize DNA object
    dna = dnaMD.dnaEY(totalBP, filename=inputFile, startBP=firstBP)

    if args.outputFile is not None:
        towrite = ''
        out = ''
        mean, result = dna.calculateLocalElasticity(bp, helical=args.helical, unit='kT')
        for i in range(result.shape[0]):
            for j in range(result.shape[0]):
                if j != result.shape[0]-1:
                    out += '{0:>10.5f}  '.format(result[i][j])
                    towrite += '{0:.5f},'.format(result[i][j])
                else:
                    out += '{0:>10.5f}\n'.format(result[i][j])
                    towrite += '{0:.5f}\n'.format(result[i][j])

        with open(args.outputFile, 'w') as fout:
            fout.write(towrite)

        print('=========== ============== Elastic Matrix =============== ===========\n')
        print(out)
        print('=========== ====================== ====================== ===========\n')

    if args.vsTime is not None:
        dna.getLocalElasticityByTime(bp, args.frameGap, helical=args.helical, unit='kT', outFile=args.vsTime)

    if args.outputFileSegments is not None:
        dna.calculateLocalElasticitySegments(bp, span=args.span, frameGap=args.frameGap, helical=args.helical,
                                             unit='kT', err_type=args.err_type, tool=args.tool,
                                             outFile=args.outputFileSegments)


def parseArguments():
    parser = argparse.ArgumentParser(prog='dnaMD localElasticity',
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

    parser.add_argument('-os', '--output-segments', action='store',
                        type=str, metavar='segments_elasticity.csv',
                        dest='outputFileSegments',
                        help=outputFileSegmentsHelp)

    parser.add_argument('-tbp', '--total-bp', action='store',
                        type=int, metavar='total-bp-number',
                        dest='totalBP', help=totalBpHelp)

    parser.add_argument('-helical', '--helical', action='store_true',
                        dest='helical', default=False,
                        help=helicalHelp)

    parser.add_argument('-bs', '--bps-start', action='store',
                        type=int, metavar='basepair-step-start-number',
                        default=1,
                        dest='startBPS', help=bpsStartHelp)

    parser.add_argument('-be', '--bps-end', action='store',
                        type=int, dest='endBPS',
                        metavar='basepair-step-end-number',
                        help=bpsEndHelp)

    parser.add_argument('-span', '--span', action='store',
                        type=int, metavar='segment-span', default=2,
                        dest='span', help=spanHelp)

    parser.add_argument('-fgap', '--frame-gap', action='store',
                        type=int, dest='frameGap',
                        metavar='frames-to-skip',
                        help=frameGapHelp)

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

    idx = sys.argv.index("localElasticity") + 1
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
