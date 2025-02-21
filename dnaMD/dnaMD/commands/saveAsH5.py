#!/usr/bin/env python
#
#
# This file is part of do_x3dna
#
# Author: Rajendra Kumar
# Copyright (C) 2014-2025  Rajendra Kumar
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
Save parameters to a HDF5 file.

The parameters can be read from the do_x3dna output files and stored in a HDF5
file. This step enables rapid reading and processing for subsequent analysis of
data.

Several calculations such as determining global helical axis, its curvature and
bending angles needs this HDF5 file because helical axis, curvature and
tangents data will be stored in this HDF5 file.

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

inputFilesHelp=\
"""List of input files containing parameters.
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

"""


outputFileHelp=\
""" Name of output file.
It is a hdf5 format file. Therefore, use ".h5" extension with it.

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

    if args.inputFiles is None:
        showErrorAndExit(parser, "No input file...\n")
    else:
        for f in args.inputFiles:
            if not os.path.isfile(f):
                showErrorAndExit(parser, "File {0} not found...\n".format(f))

    # Determine file-extension type
    outputFileExtension = os.path.splitext(args.outputFile)[1]
    if outputFileExtension not in ['.h5', '.hdf5', '.hdf']:
        showErrorAndExit(parser, "File extension {0} is not recognized as an \
        acceptable HDF5 extension.\n Use '.h5', '.hdf5' \
        or '.hdf'.".format(outputFileExtension))

    # initialize DNA object
    dna = dnaMD.DNA(totalBP, filename=args.outputFile, startBP=firstBP)
    for f in args.inputFiles:
        dnaMD.setParametersFromFile(dna, f)
    del dna

def parseArguments():
    parser = argparse.ArgumentParser(prog='dnaMD saveAsH5',
                                    description=description,
                                    formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-tbp', '--total-bp', action='store',
                        type=int, metavar='total-bp-number',
                        dest='totalBP', help=totalBpHelp)

    parser.add_argument('-fbp', '--first-bp', action='store',
                        type=int, metavar='1', default=1,
                        dest='firstBP', help=bpFirstHelp)

    parser.add_argument('-i', '--input-files', action='store',
                        type=lambda s: [item for item in s.split(',')],
                        metavar='input1.dat,input2.dat,input3.dat',
                        dest='inputFiles', help=inputFilesHelp)

    parser.add_argument('-o', '--output', action='store',
                        type=str, metavar='output.h5',
                        dest='outputFile', help=outputFileHelp)

    idx = sys.argv.index("saveAsH5")+1
    args = parser.parse_args(args=sys.argv[idx:])

    return parser, args

def showErrorAndExit(parser, message):
    parser.print_help()
    print("\n===== ERROR =======")
    print(message)
    print("See Usage Above!!!")
    sys.exit(False)
