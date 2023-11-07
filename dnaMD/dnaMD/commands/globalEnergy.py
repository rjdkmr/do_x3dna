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
import argparse
import numpy as np

import dnaMD

description="""DESCRIPTION
===========
Global Deformation Energy of the DNA

This can be used to Global Deformation Energy of the DNA from the simulations. At first, elastic matrix from reference
DNA (most often free or unbound DNA) is calculated and subsequently this matrix is used to calculate deformation free 
energy of probe DNA (most often bound DNA).

"""

inputRefFileHelp=\
"""Name of input reference file (hdf5 file).
File containing parameters of reference DNA for which global elastic properties will
be calculated. Most often it is free or unbound DNA.

This file should contain the required parameters. It should be hdf5 storage file.

"""

inputProbeFileHelp=\
"""Name of input probe file (hdf5 file).
File containing parameters of probe DNA for which global deformation energy will
be calculated. Most often it is bound DNA.

This file should contain the required parameters. It should be hdf5 storage file.

"""

outputFileHelp=\
"""Name of output file in csv format.
This file will contain the energy values as a function of time.

"""

energyTermHelp=\
"""Energy terms to be calculated.
For which motions, energy should be calculated.

Following keywords are available:
* all         : (Default) All below listed energy terms will be calculated
* full        : Use entire elastic matrix -- all motions with their coupling
* diag        : Use diagonal of elastic matrix -- all motions but no coupling
* b1          : Only bending-1 motion
* b2          : Only bending-2 motion
* stretch     : Only stretching motion
* twist       : Only Twisting motions
* st_coupling : Only stretch-twist coupling motion
* bs_coupling : Only Bending and stretching coupling
* bt_coupling : Only Bending and Twisting coupling
* bb_coupling : Only bending-1 and bending-2 coupling
* bend        : Both bending motions with their coupling
* st          : Stretching and twisting motions with their coupling
* bs          : Bending (b1, b2) and stretching motions with their coupling
* bt          : Bending (b1, b2) and twisting motions with their coupling

In case of elasticity type "ST", only following four energy terms are available "all", "diag", "stretch", "twist" and
"st_coupling".

The terms should provided as comma separated values. e.g. -et "full,diag,b1,b2,stretch,twist".

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
" "ST"  : Stretching-Twisting --- Bending motions are ignored.

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

enGlobalTypes = ['full', 'diag',  'stretch', 'twist', 'st_coupling', 'b1', 'b2',
                      'bend', 'bs_coupling', 'bt_coupling', 'bb_coupling', 'st', 'bs', 'bt' ]


def main():
    parser, args = parseArguments()

    # Input file
    inputRefFile = None
    if args.inputRefFile is not None:
        inputRefFile = args.inputRefFile.name
        args.inputRefFile.close()
    else:
        showErrorAndExit(parser, "No Input File for Reference DNA!!!\n")

    # Input file
    inputProbeFile = None
    if args.inputProbeFile is not None:
        inputProbeFile = args.inputProbeFile.name
        args.inputProbeFile.close()
    else:
        showErrorAndExit(parser, "No Input File for Probe DNA!!!\n")

    # Determine file-extension type
    inputFileExtension = os.path.splitext(inputRefFile)[1]
    if inputFileExtension not in ['.h5', '.hdf5', '.hdf']:
        showErrorAndExit(parser, "Input file for Reference DNA should be in HDF5 (h5, hdf5 or hdf extension) format.\n")
    inputFileExtension = os.path.splitext(inputProbeFile)[1]
    if inputFileExtension not in ['.h5', '.hdf5', '.hdf']:
        showErrorAndExit(parser, "Input file for probe DNA should be in HDF5 (h5, hdf5 or hdf extension) format.\n")

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

    # Check energy terms and make a list
    outEnergyTerms = checkEnergyTerms(args)

    # initialize DNA object
    dna = dnaMD.dnaEY(totalBP, esType=args.esType, filename=inputRefFile, startBP=firstBP)
    complexDna = dnaMD.DNA(totalBP, filename=inputProbeFile, startBP=firstBP)

    # Check if mask is in object
    if dna.dna.mask is not None:
        masked = True
    else:
        masked = False

    time, energy = dna.getGlobalDeformationEnergy(bp, complexDna, paxis=args.paxis, which=outEnergyTerms, masked=masked,
                                                  outFile = args.outputFile)

    if args.err_type is not None:
        error = dnaMD.get_error(time, list(energy.values()), len(outEnergyTerms), err_type=args.err_type, tool=args.tool)

        sys.stdout.write("==============================================\n")
        sys.stdout.write('{0:<16}{1:>14}{2:>14}\n'.format('Energy(kJ/mol)', 'Average', 'Error'))
        sys.stdout.write("----------------------------------------------\n")
        for i in range(len(outEnergyTerms)):
            sys.stdout.write('{0:<16}{1:>14.3f}{2:>14.3f}\n'.format(outEnergyTerms[i], np.mean(energy[outEnergyTerms[i]]),
                                                                    error[i]))
        sys.stdout.write("==============================================\n\n")


def checkEnergyTerms(args):
    if args.esType == 'BST':
        energyTerms = enGlobalTypes
    else:
        energyTerms = enGlobalTypes[:5]

    outEnergyTerms = args.energyTerms
    if 'all' in outEnergyTerms:
        outEnergyTerms = energyTerms
    else:
        for key in outEnergyTerms:
            if key not in energyTerms:
                raise ValueError('{0} is not a supported keyword.\n Use from the following list: \n{1}'.format(
                    outEnergyTerms, energyTerms))

    return outEnergyTerms

def parseArguments():
    parser = argparse.ArgumentParser(prog='dnaMD globalEnergy',
                                     description=description,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-ir', '--input-ref', action='store',
                        type=argparse.FileType('rb'), metavar='ref_dna.h5',
                        dest='inputRefFile', required=False, help=inputRefFileHelp)

    parser.add_argument('-ip', '--input-probe', action='store',
                        type=argparse.FileType('rb'), metavar='probe_dna.h5',
                        dest='inputProbeFile', required=False, help=inputProbeFileHelp)

    parser.add_argument('-o', '--output', action='store',
                        type=str, metavar='output.dat',
                        dest='outputFile', help=outputFileHelp)

    parser.add_argument('-et', '--energy-terms', action='store',
                        type=lambda s: [item.rstrip().lstrip() for item in s.split(',')],
                        metavar='"full,diag,strecth,twist"', default='all',
                        dest='energyTerms', help=energyTermHelp)

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

    parser.add_argument('-paxis', '--principal-axis', action='store',
                        type=str, metavar='X', default=None,
                        choices=['X', 'Y', 'Z'],
                        dest='paxis', help=paxisHelp)

    parser.add_argument('-em', '--error-method', action='store',
                        type=str, metavar='block', default='block',
                        choices=['std', 'acf', 'block'],
                        dest='err_type', help=errorHelp)

    parser.add_argument('-gt', '--gromacs-tool', action='store',
                        type=str, metavar='gmx analyze', default='gmx analyze',
                        dest='tool', help=toolsHelp)

    parser.add_argument('-fbp', '--first-bp', action='store',
                        type=int, metavar='1', default=1,
                        dest='firstBP', help=bpFirstHelp)

    idx = sys.argv.index("globalEnergy") + 1
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
