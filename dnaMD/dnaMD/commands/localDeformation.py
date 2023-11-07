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

description= """DESCRIPTION
===========
Deformation of local parameters in probe DNA with respect to a reference DNA

This can be used to calculate deformation in the given parameters of a probe DNA with respect to a reference DNA along 
the base-pairs/steps.

NOTE: Number of input segments/bp/bps should match between probe and reference DNA. It means ([-bePrb] - [-bsPrb])
should be equal to ([-beRef] - [-bsRef]).

"""

probeInputFileHelp=\
"""Name of input file for probe DNA (from do_x3dna or hdf5 file).
This file should contain the required parameters. It can be a file either
produced from do_x3dna or hdf5 storage file. 

"""

refInputFileHelp=\
"""Name of input file for reference DNA (from do_x3dna or hdf5 file).
This file should contain the required parameters. It can be a file either
produced from do_x3dna or hdf5 storage file. 

"""

outputFileHelp=\
"""Name of output file.
The extracted output will be written in output file.

"""

xaxisBPHelp=\
""" Basepair/step number on the x-axis in output file
To write the base-pair/step number of either probe or reference DNA in the output file.

If base-pair numbering is different in reference and probe DNAs, then this option can be used to
write the base-pair number of either probe (-oxb probe) or reference (-oxb ref) DNA in the output file.

"""

parameterHelp=\
"""Parameter name.
This parameter will be extracted from file. Ensure that parameter is present
in the file, otherwise wrong values will be extracted from file.

"""

probeTotalBPHelp=\
"""Total number of basepair in probe DNA.
It is an essential input.

"""

refTotalBPHelp=\
"""Total number of basepair in reference DNA.
It is an essential input.

"""

probeFirstBPHelp=\
"""Basepair number of first base-pair in probe DNA.
Usually it is one. Therefore, if this option is not provided, base-pair
numbering will start from one.

"""

refFirstBPHelp=\
"""Basepair number of first base-pair in reference DNA.
Usually it is one. Therefore, if this option is not provided, base-pair
numbering will start from one.

"""

probeStartBPHelp=\
"""First BP/BPS of probe DNA after which parameter will be extracted.
If it is not given, first basepair or base-step will be considered.

"""

probeEndBPHelp=\
"""Last BP/BPS of probe DNA upto which parameter will be extracted.
If it is not given, last basepair or base-step will be considered.

"""

refStartBPHelp=\
"""First BP/BPS of reference DNA after which parameter will be extracted.
If it is not given, first basepair or base-step will be considered.

"""

refEndBPHelp=\
"""Last BP/BPS of reference DNA upto which parameter will be extracted.
If it is not given, last basepair or base-step will be considered.

"""

mergeMethodHelp=\
"""Method to merge the parameters of consecutive basepairs/steps given by -mb/--merge-bps.

Currently accepted keywords are as follows:
    * mean : Average of local parameters
    * sum : Sum of local parameters

"""

errorHelp=\
""" Error method
It can be as following:
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

merge_bp_help=\
"""Number of consecutive base-pairs/steps to merge for creating the small non-overlapping DNA
segments. By default, averages and errors are calculated for each base-pair/step separately.
However, averages  and errors can also  be calculated for small non-overlapping DNA segment
by merging the parameters of consecutive base-pairs/steps.

"""


def main():
    parser, args = parseArguments()

    # Input file
    probeInputFile = None
    if args.probeInputFile is not None:
        probeInputFile = args.probeInputFile.name
        args.probeInputFile.close()
    else:
        showErrorAndExit(parser, "No Probe DNA Input File!!!\n")
    refInputFile = None
    if args.refInputFile is not None:
        refInputFile = args.refInputFile.name
        args.refInputFile.close()
    else:
        showErrorAndExit(parser, "No reference DNA Input File!!!\n")

    # Determine file-extension type
    probeFileType = 'text'
    inputFileExtension = os.path.splitext(probeInputFile)[1]
    if inputFileExtension in ['.h5', '.hdf5', 'hdf']:
        probeFileType = 'hdf5'

    refFileType = 'text'
    inputFileExtension = os.path.splitext(refInputFile)[1]
    if inputFileExtension in ['.h5', '.hdf5', 'hdf']:
        refFileType = 'hdf5'

    # Total number of base-pair  --- probe
    probeFirstBP = args.probeFirstBP
    probeTotalBP = None
    if args.probeTotalBP is None:
        showErrorAndExit(parser, "No total number of BP for probe DNA!!!\n")
    else:
        probeTotalBP = args.probeTotalBP

    # Total number of base-pair  --- reference
    refFirstBP = args.probeFirstBP
    refTotalBP = None
    if args.refFirstBP is None:
        showErrorAndExit(parser, "No total number of BP for reference DNA!!!\n")
    else:
        refTotalBP = args.refTotalBP

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

    # Determine start and end-bp ---- probe
    probeStartBP = args.probeStartBP
    if probeStartBP is None:
        probeStartBP = probeFirstBP

    probeEndBP = None
    if args.probeEndBP is not None:
        probeEndBP = args.probeEndBP
    else:
        probeEndBP = probeTotalBP - toMinusBP

    probeBP = [probeStartBP, probeEndBP]

    # Determine start and end-bp ---- reference
    refStartBP = args.refStartBP
    if refStartBP is None:
        refStartBP = refFirstBP

    refEndBP = None
    if args.refEndBP is not None:
        refEndBP = args.refEndBP
    else:
        refEndBP = refTotalBP - toMinusBP

    refBP = [refStartBP, refEndBP]

    # Check requested segment is equal in length
    if (refEndBP- refStartBP) != (probeEndBP - probeStartBP):
        msg='Length of requested segment is not equal in probe ({0}) and reference ({1}) DNA.'.format(probeBP, refBP)
        showErrorAndExit(parser, msg)


    # Check consistency of start bp --- probe
    if (probeStartBP < probeFirstBP) or (probeStartBP > probeTotalBP+probeFirstBP-toMinusBP):
        msg = 'The requested start bp {0} is out side of {1}-{2} range.'.format(probeStartBP, probeFirstBP, probeTotalBP+probeFirstBP-toMinusBP)
        showErrorAndExit(parser, msg)

    # Check consistency of start bp --- reference
    if (refStartBP < refFirstBP) or (refStartBP > refTotalBP+refFirstBP-toMinusBP):
        msg = 'The requested start bp {0} is out side of {1}-{2} range.'.format(refStartBP, refFirstBP, refTotalBP+refFirstBP-toMinusBP)
        showErrorAndExit(parser, msg)

    # Check consistency of end-bp  --- probe
    if probeStartBP > probeEndBP:
        msg = 'The requested end bp {0} is larger than requested start bp {1}!!!'.format(probeEndBP, probeStartBP)
        showErrorAndExit(parser, msg)
    if (probeEndBP > probeTotalBP+probeFirstBP-toMinusBP) or (probeEndBP < probeFirstBP):
        msg = 'The requested end bp {0} is out side of {1}-{2} range.'.format(probeEndBP, probeFirstBP, probeTotalBP+probeFirstBP-toMinusBP)
        showErrorAndExit(parser, msg)

    # Check consistency of end-bp  --- reference
    if refStartBP > refEndBP:
        msg = 'The requested end bp {0} is larger than requested start bp {1}!!!'.format(refEndBP, refStartBP)
        showErrorAndExit(parser, msg)
    if (refEndBP > refTotalBP+refFirstBP-toMinusBP) or (refEndBP < refFirstBP):
        msg = 'The requested end bp {0} is out side of {1}-{2} range.'.format(refEndBP, refFirstBP, refTotalBP+refFirstBP-toMinusBP)
        showErrorAndExit(parser, msg)

    # Check if merge-method is given
    if args.merge_method is None:
        msg = 'No merging method is provided!!!!'
        showErrorAndExit(parser, msg)

    # Store input file name
    probeFilename = None
    if probeFileType == 'hdf5':
        probeFilename = probeInputFile
    refFilename = None
    if refFileType == 'hdf5':
        refFilename = refInputFile

    # Check merge_bp
    if args.merge_bp > probeTotalBP or args.merge_bp > refTotalBP:
        msg = ' Number of bp/s to merge is larger than total number of bp.'
        showErrorAndExit(parser, msg)

    # check gromacs tools
    if 'analyze' not in args.tool:
        msg = '{0} might not be suitable for error calculation. \n Use Gromacs analyze tool g_analyze or "gmx analyze".'\
            .format(args.tool)
        showErrorAndExit(parser, msg)

    # initialize DNA object
    dnaProbe = dnaMD.DNA(probeTotalBP, filename=probeInputFile, startBP=probeFirstBP)
    dnaRef = dnaMD.DNA(refTotalBP, filename=refInputFile, startBP=refFirstBP)

    # Check if mask is in object
    if (dnaProbe.mask is not None) and (args.parameter in dnaMD.maskedParameters):
        masked = True
    else:
        masked = False

    # Directly read the data from text file  --- probe
    if probeFileType == 'text':
        # No support for smoothed axis, curvature and tangent
        if args.parameter in ['helical x-axis smooth', 'helical y-axis smooth', 'helical z-axis smooth', 'helical axis curvature', 'helical axis tangent']:
            print("\n===== ERROR =======")
            print('Extraction of parameter "{0}" is only supported via input HDF5 file.'.format(args.parameter))
            print("\n===== ERROR =======")
            sys.exit(-1)
        # Read and load parameter from file
        dnaMD.setParametersFromFile(dnaProbe, probeInputFile, args.parameter, bp=probeBP)

    # Directly read the data from text file  --- reference
    if refFileType == 'text':
        # No support for smoothed axis, curvature and tangent
        if args.parameter in ['helical x-axis smooth', 'helical y-axis smooth', 'helical z-axis smooth', 'helical axis curvature', 'helical axis tangent']:
            print("\n===== ERROR =======")
            print('Extraction of parameter "{0}" is only supported via input HDF5 file.'.format(args.parameter))
            print("\n===== ERROR =======")
            sys.exit(-1)
        # Read and load parameter from file
        dnaMD.setParametersFromFile(dnaRef, refInputFile, args.parameter, bp=refBP)

    # Main calculation here
    bpOutRef, bpOutProbe, deformation, error = dnaMD.localDeformationVsBPS(dnaRef, refBP, dnaProbe, probeBP,
                                                                           args.parameter,
                                                                           err_type=args.err_type,
                                                                           bp_range=True, merge_bp=args.merge_bp,
                                                                           merge_method=args.merge_method,
                                                                           masked=masked, tool=args.tool)

    # Write the extracted data in a text file
    fout = open(args.outputFile, 'w')
    fout.write('# bp(mid) \t {0}-avg \t {0}-error\n'.format(args.parameter))
    for i in range(len(deformation)):
        if args.xaxisBP == 'probe':
            fout.write("{0}\t".format(bpOutProbe[i]))
        else:
            fout.write("{0}\t\t".format(bpOutRef[i]))
        fout.write("{0:.6}\t{1:.6}\n".format(deformation[i], error[i]))
    fout.close()


def parseArguments():
    parser = argparse.ArgumentParser(prog='dnaMD localDeformation',
                                    description=description,
                                    formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-ir', '--input-reference', action='store',
                        type=argparse.FileType('rb'), metavar='L-BP_fdna.dat',
                        dest='refInputFile', required=False, help=refInputFileHelp)

    parser.add_argument('-ip', '--input-probe', action='store',
                        type=argparse.FileType('rb'), metavar='L-BP_cdna.dat',
                        dest='probeInputFile', required=False, help=probeInputFileHelp)

    parser.add_argument('-o', '--output', action='store',
                        type=str, metavar='output.dat',
                        dest='outputFile', help=outputFileHelp)

    parser.add_argument('-oxb', '--xaxis-bp-source', action='store',
                        type=str, metavar='ref', default='ref',
                        choices=['probe', 'ref'],
                        dest='xaxisBP', help=xaxisBPHelp)

    parser.add_argument('-tbpPrb', '--total-bp-probe', action='store',
                        type=int, metavar='total-bp-number',
                        dest='probeTotalBP', help=probeTotalBPHelp)

    parser.add_argument('-tbpRef', '--total-bp-Ref', action='store',
                        type=int, metavar='total-bp-number',
                        dest='refTotalBP', help=refTotalBPHelp)

    parser.add_argument('-p', '--parameter', action='store',
                        dest='parameter', metavar='parameter',
                        type=str, help=parameterHelp)

    parser.add_argument('-em', '--error-method', action='store',
                        type=str, metavar='block', default='block',
                        choices=['std', 'acf', 'block'],
                        dest='err_type', help=errorHelp)

    parser.add_argument('-mb', '--merge-bps', action='store',
                        type=int, metavar='bp/s number', default=1,
                        dest='merge_bp', help=merge_bp_help)

    parser.add_argument('-gt', '--gromacs-tool', action='store',
                        type=str, metavar='g_analyze', default='g_analyze',
                        dest='tool', help=toolsHelp)

    parser.add_argument('-bsPrb', '--bp-start-probe', action='store',
                        type=int, metavar='bp/s-start-number',
                        dest='probeStartBP', help=probeStartBPHelp)

    parser.add_argument('-bePrb', '--bp-end-probe', action='store',
                        type=int, dest='probeEndBP',
                        metavar='bp/s-end-number',
                        help=probeEndBPHelp)

    parser.add_argument('-bsRef', '--bp-start-ref', action='store',
                        type=int, metavar='bp/s-start-number',
                        dest='refStartBP', help=refStartBPHelp)

    parser.add_argument('-beRef', '--bp-end-ref', action='store',
                        type=int, dest='refEndBP',
                        metavar='bp/s-end-number',
                        help=refEndBPHelp)

    parser.add_argument('-mm', '--merge-method', action='store',
                        type=str, dest='merge_method',
                        metavar='sum-or-mean',
                        choices=['mean', 'sum'],
                        help=mergeMethodHelp)

    parser.add_argument('-fbpPrb', '--first-bp-probe', action='store',
                        type=int, metavar='1', default=1,
                        dest='probeFirstBP', help=probeFirstBPHelp)

    parser.add_argument('-fbpRef', '--first-bp-ref', action='store',
                        type=int, metavar='1', default=1,
                        dest='refFirstBP', help=refFirstBPHelp)

    idx = sys.argv.index("localDeformation")+1
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
