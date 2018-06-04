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
Parameter distribution during simulation

This can be used to calculate distribution of a parameter of either a specific base-pair/step
or over a DNA segment during the simulations.

"""

inputFileHelp=\
"""Name of input file (from do_x3dna or hdf5 file).
This file should contain the required parameters. It can be a file either
produced from do_x3dna or hdf5 storage file.

"""

outputFileHelp=\
"""Name of output file.
The extracted output will be written in output file.

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
parameterHelp=\
"""Parameter name.
This parameter will be extracted from file. Ensure that parameter is present
in the file, otherwise wrong values will be extracted from file.

"""

bpStartHelp=\
"""First BP/BPS of DNA after which parameter will be extracted.
If it is not given, first basepair or base-step will be considered.

"""

bpEndHelp=\
"""Last BP/BPS of DNA upto which parameter will be extracted.

If it is not given, parameter for only a single bp/s given with -bs/--bp-start
option will be extracted.

"""

mergeMethodHelp=\
"""Method to merge the parameter of a DNA segment from local parameters
of all base-pairs/steps that are within the range given by '-bs' and '-be'.

Currently accepted keywords are as follows:
    * mean : Average of local parameters
    * sum : Sum of local parameters

When only "-bs" option is provided without "-be", then -mm/--merge-method is
not required.

"""

binsHelp=\
""" Number of bins in the histogram
Default value is 30.

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
    fileType = 'text'
    inputFileExtension = os.path.splitext(inputFile)[1]
    if inputFileExtension in ['.h5', '.hdf5', 'hdf']:
        fileType = 'hdf5'

    # Total number of base-pair
    firstBP = args.firstBP
    totalBP = None
    if args.totalBP is None:
        showErrorAndExit(parser, "No total number of BP!!!\n")
    else:
        totalBP = args.totalBP

    # Check for input parameter
    if args.parameter is None:
        showErrorAndExit(parser, "No Parameter name!!!\n")

    if not ( (args.parameter in dnaMD.basePairTypeParameters) or \
            (args.parameter in dnaMD.baseStepTypeParameters) ):
        parser.print_help()
        print("\n===== ERROR =======")
        print('Unknown parameter name "{0}"!!!\n'.format(args.parameter))
        print("Accepted parameters are as follows:")
        count = 1
        for parameter in dnaMD.basePairTypeParameters:
            print('{0}. "{1}"'.format(count, parameter))
            count += 1
        for parameter in dnaMD.baseStepTypeParameters:
            print('{0}. "{1}"'.format(count, parameter))
            count += 1
        print("\n===================")
        sys.exit(-1)

    # Output file
    if args.outputFile is None:
        showErrorAndExit(parser, "No output File!!!\n")

    paramType = dnaMD.getParameterType(args.parameter)
    toMinusBP = 1
    if paramType == 'bps':
        toMinusBP = 2

    # Determine start and end-bp
    startBP = args.startBP
    if startBP is None:
        startBP = firstBP
    endBP = None
    if args.endBP is not None:
        endBP = args.endBP

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

    # Check if merge is required
    if endBP is None:
        merge = False
    else:
        merge = True

    # Check if merge-method is given
    if merge and args.merge_method is None:
        msg = 'No merging method is provided!!!!'
        showErrorAndExit(parser, msg)

    # Store input file name
    filename = None
    if fileType == 'hdf5':
        filename = inputFile

    if endBP is not None:
        bp = [startBP, endBP]
    else:
        bp = [startBP]

    # initialize DNA object
    dna = dnaMD.DNA(totalBP, filename=filename, startBP=firstBP)

    # Check if mask is in object
    if dna.mask is not None:
        masked = True
    else:
        masked = False

    # Directly read the data from text file
    if fileType == 'text':
        # No support for smoothed axis, curvature and tangent
        if args.parameter in ['helical x-axis smooth', 'helical y-axis smooth', 'helical z-axis smooth', 'helical axis curvature', 'helical axis tangent']:
            print("\n===== ERROR =======")
            print('Extraction of parameter "{0}" is only supported via input HDF5 file.'.format(args.parameter))
            print("\n===== ERROR =======")
            sys.exit(-1)

        # Read and load parameter from file
        dnaMD.setParametersFromFile(dna, inputFile, args.parameter, bp=bp)

    # Extract the input parameter for input DNA/RNA segment
    values, density = dna.parameter_distribution(args.parameter, bp, bins=30, merge=merge, merge_method=args.merge_method, masked=masked)

    # Write the extracted data in a text file
    fout = open(args.outputFile, 'w')
    fout.write('# "{0}" \t Density\n'.format(args.parameter))
    for i in range(len(values)):
        fout.write("{0:.6}\t{1:.6}\n".format(values[i], density[i]))
    fout.close()

def parseArguments():
    parser = argparse.ArgumentParser(prog='dnaMD histogram',
                                    description=description,
                                    formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', action='store',
                        type=argparse.FileType('rb'), metavar='L-BP_cdna.dat',
                        dest='inputFile', required=False, help=inputFileHelp)

    parser.add_argument('-o', '--output', action='store',
                        type=str, metavar='output.dat',
                        dest='outputFile', help=outputFileHelp)

    parser.add_argument('-tbp', '--total-bp', action='store',
                        type=int, metavar='total-bp-number',
                        dest='totalBP', help=totalBpHelp)

    parser.add_argument('-bins', '--bins', action='store',
                        type=int, metavar='bins', default=30,
                        dest='bins', help=binsHelp);

    parser.add_argument('-p', '--parameter', action='store',
                        dest='parameter', metavar='parameter',
                        type=str, help=parameterHelp)

    parser.add_argument('-bs', '--bp-start', action='store',
                        type=int, metavar='bp/s-start-number',
                        default=1,
                        dest='startBP', help=bpStartHelp)

    parser.add_argument('-be', '--bp-end', action='store',
                        type=int, dest='endBP',
                        metavar='bp/s-end-number',
                        help=bpEndHelp)

    parser.add_argument('-mm', '--merge-method', action='store',
                        type=str, dest='merge_method',
                        metavar='sum-or-mean',
                        choices=['mean', 'sum'],
                        help=mergeMethodHelp)

    parser.add_argument('-fbp', '--first-bp', action='store',
                        type=int, metavar='1', default=1,
                        dest='firstBP', help=bpFirstHelp)

    idx = sys.argv.index("histogram") + 1
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
