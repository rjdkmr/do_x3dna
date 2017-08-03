import os
import sys
import re
import argparse

import dnaMD

description=\
"""Paramter as a function of time

This can be used to extract the parameter of either a specfic base-pair/step
or over a DNA segment as a function of time.

"""

inputFileHelp=\
"""Name of input file.
This file should contain the required parameters.

"""

ouputFileHelp=\
"""Name of ouput file.
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
This paramter will be extracted from file. Ensure that parameter is present
in the file, otherwise wrong values will be extracted from file.

"""

bpStartHelp=\
"""First BP/BPS of a segment for extraction
If it is not given, first basepair or base-step will be considered.

"""

bpEndHelp=\
"""Last BP/BPS of a segment for extraction

If it is not given, parameter for only a single bp/s given with -bs/--bp-start
option will be extracted.

"""

mergeMethodHelp=\
"""Method to merge the paremeter of a DNA segment from local parameters
of all base-pairs/steps that are within the range given by '-bs' and '-be'.

Currently accepted keywords are as follows:
    * mean : Average of local parameters
    * sum : Sum of local parameters

When only "-bs" option is provided without "-be", then -mm/--merge-method is
not required.

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
    inputFileExtension = os.path.splitext(inputFile)
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

    if not ( (args.parameter in dnaMD.basepairParameterNames) or \
            (args.parameter in dnaMD.basestepParamterNames) ):
        parser.print_help()
        print("\n===== ERROR =======")
        print('Unknown parameter name "{0}"!!!\n'.format(args.parameter))
        print("Accepted parameters are as follows:")
        count = 1
        for parameter in dnaMD.basepairParameterNames:
            print('{0}. "{1}"'.format(count, parameter))
            count += 1
        for parameter in dnaMD.basestepParamterNames:
            print('{0}. "{1}"'.format(count, parameter))
            count += 1
        print("\n===================")
        sys.exit(-1)

    # Output file
    if args.ouputFile is None:
        showErrorAndExit(parser, "No output File!!!\n")

    paramType = dnaMD.getParameterType(args.parameter)
    toMinusBP = 1
    if paramType == 'bps':
        toMinusBP = 2

    # Determine start and end-bp
    startBP = args.startBP
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
    if dna.mask:
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
    time, value = dna.time_vs_parameter(args.parameter, bp, merge=merge, merge_method=args.merge_method, masked=masked)

    # Write the extracted data in a text file
    fout = open(args.ouputFile, 'w')
    fout.write('# Time \t "{0}"\n'.format(args.parameter))
    for i in range(len(time)):
        fout.write("{0}\t{1}\n".format(time[i], value[i]))
    fout.close()

def parseArguments():
    parser = argparse.ArgumentParser(prog='dnaMD vsTime',
                                    description=description,
                                    formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', action='store',
                        type=argparse.FileType('rb'), metavar='L-BP_cdna.dat',
                        dest='inputFile', required=False, help=inputFileHelp)

    parser.add_argument('-o', '--output', action='store',
                        type=str, metavar='output.dat',
                        dest='ouputFile', help=ouputFileHelp)

    parser.add_argument('-tbp', '--total-bp', action='store',
                        type=int, metavar='total-bp-number',
                        dest='totalBP', help=totalBpHelp)

    parser.add_argument('-fbp', '--first-bp', action='store',
                        type=int, metavar='1', default=1,
                        dest='firstBP', help=bpFirstHelp)

    parser.add_argument('-p', '--parameter', action='store',
                        dest='parameter', metavar='parameter',
                        type=str, help=parameterHelp)

    parser.add_argument('-bs', '--bp-start', action='store',
                        type=int, metavar='bp/s-start-number',
                        default=1,
                        dest='startBP', help=bpStartHelp)

    parser.add_argument('-be', '--bp-end', action='store',
                        type=int, dest='endBP',
                        metavar='b/ps-end-number',
                        help=bpEndHelp)

    parser.add_argument('-mm', '--merge-method', action='store',
                        type=str, dest='merge_method',
                        metavar='sum-or-mean',
                        choices=['mean', 'sum'],
                        help=mergeMethodHelp)

    idx = sys.argv.index("vsTime")+1
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
