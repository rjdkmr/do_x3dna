import os
import sys
import re
import argparse

import dnaMD

description=\
"""Save parameters to a HDF5 file
=================================
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

lbpInFileHelp = \
"""File containing base-pair parameters.
This file is obtained from do_x3dna. It should contain six parameters:
(1) shear, (2) stretch, (3) stagger, (4) buckle, (5) propeller and (6) opening.

"""

lbpsInFileHelp=\
""" File containing base-step parameters.
This file is obtained from do_x3dna. It should contain six parameters:
(1) shift, (2) slide, (3) rise, (4) tilt, (5) roll and (6) twist

"""

lbphInFileHelp=\
""" File containing helical base-step parameters.
This file is obtained from do_x3dna. It should contain six parameters:
(1) x-disp, (2) y-disp, (3) h-rise, (4) inclination, (5) tip and (6) h-twist

"""

bbDihInFileHelp=\
""" File containing backbone and chi dihedral angles.
This file is obtained from do_x3dna. It should contain 14 parameters:
(1) alpha S1, (2) beta S1, (3) gamma S1, (4) delta S1, (5) epsilon S1,
(6) zeta S1,  (7) chi S1,  (8) alpha S2, (9) beta S2,  (10) gamma S2,
(11) delta S2, (12) epsilon S2, (13)zeta S2, and (14) chi S2

where S1 and S2 is first and second strand respectively.

"""

groovesInFileHelp=\
""" File containing major and minor grooves width.
This file is obtained from do_x3dna. It should contain four parameters:
(1) minor groove, (2) minor groove refined,
(3) major groove, (4) major groove refined.

"""

haxisInFileHelp=\
""" File containing local helical axis coordinates.
This file is obtained from do_x3dna. It should contain coordinates of local
helical axis.

Note that global and smooth helical axis can be calculated using same HDF5
file as input.

"""
hradInFileHelp=\
"""File containing helical radius.
This file is obtained from do_x3dna. It should contain local helical radius.

"""
ouputFileHelp=\
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

    inputFilesDict = checkForInputFile(args, parser)

    # Determine file-extension type
    fileType = 'hdf5'
    ouputFileExtension = os.path.splitext(args.ouputFile)[1]
    if ouputFileExtension not in ['.h5', '.hdf5', 'hdf']:
        showErrorAndExit(parser, "File extension {0} is not recognized as an \
        acceptable HDF5 extension.\n Use '.h5', '.hdf5' \
        or '.hdf'.".format(ouputFileExtension))

    # initialize DNA object
    dna = dnaMD.DNA(totalBP, filename=args.ouputFile, startBP=firstBP)

    if 'bp' in inputFilesDict:
        dnaMD.setParametersFromFile(dna, inputFilesDict['bp'], dnaMD.basePairParameters[:])
    if 'bps' in inputFilesDict:
        dnaMD.setParametersFromFile(dna, inputFilesDict['bps'], dnaMD.baseStepParameters[:])
    if 'bph' in inputFilesDict:
        dnaMD.setParametersFromFile(dna, inputFilesDict['bph'], dnaMD.helicalBaseStepParameters[:])
    if 'bbdih' in inputFilesDict:
        dnaMD.setParametersFromFile(dna, inputFilesDict['bbdih'], dnaMD.backboneDihedrals[:])
    if 'grooves' in inputFilesDict:
        dnaMD.setParametersFromFile(dna, inputFilesDict['grooves'], dnaMD.groovesParameters[:])
    if 'haxis' in inputFilesDict:
        dnaMD.setParametersFromFile(dna, inputFilesDict['haxis'], dnaMD.helicalAxisParameters[:])
    if 'hrad' in inputFilesDict:
        dnaMD.setParametersFromFile(dna, inputFilesDict['hrad'], dnaMD.helicalRadiusParameters[:])


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

    parser.add_argument('-lbp', '--local-base-pair', action='store',
                        type=argparse.FileType('r'), metavar='L-BP_g.dat',
                        dest='lbpInFile', help=lbpInFileHelp)

    parser.add_argument('-lbps', '--local-base-step', action='store',
                        type=argparse.FileType('r'), metavar='L-BPS_g.dat',
                        dest='lbpsInFile', help=lbpsInFileHelp)

    parser.add_argument('-lbph', '--local-base-helical-step', action='store',
                        type=argparse.FileType('r'), metavar='L-BPH_g.dat',
                        dest='lbphInFile', help=lbphInFileHelp)

    parser.add_argument('-bbd', '--backbone-dihedrals', action='store',
                        type=argparse.FileType('r'), metavar='BackBoneCHiDihedrals_g.dat',
                        dest='bbDihInFile', help=bbDihInFileHelp)

    parser.add_argument('-mg', '--major-minor-grooves', action='store',
                        type=argparse.FileType('r'), metavar='MGroove_g.dat',
                        dest='groovesInFile', help=groovesInFileHelp)

    parser.add_argument('-ha', '--helical-axis', action='store',
                        type=argparse.FileType('r'), metavar='HelAxis_g.dat',
                        dest='haxisInFile', help=haxisInFileHelp)

    parser.add_argument('-hr', '--helical-radius', action='store',
                        type=argparse.FileType('r'), metavar='HelixRad_g.dat',
                        dest='hradInFile', help=hradInFileHelp)

    parser.add_argument('-o', '--output', action='store',
                        type=str, metavar='output.h5',
                        dest='ouputFile', help=ouputFileHelp)

    idx = sys.argv.index("saveAsH5")+1
    args = parser.parse_args(args=sys.argv[idx:])

    return parser, args

def checkForInputFile(args, parser):

    allInputFilesDict = { 'bp': args.lbpInFile,
                          'bps' : args.lbpsInFile,
                          'bph' : args.lbphInFile,
                          'bbdih' : args.bbDihInFile,
                          'grooves' : args.groovesInFile,
                          'haxis' : args.haxisInFile,
                          'hrad' : args.hradInFile }

    inputFilesDict = dict()
    for f in allInputFilesDict:
        if allInputFilesDict[f] is not None:
            inputFilesDict[f] = allInputFilesDict[f].name
            allInputFilesDict[f].close()

    if not inputFilesDict:
        msg = 'No input file provided. At least one input file is required for processing.'
        showErrorAndExit(parser, msg)

    return inputFilesDict

def showErrorAndExit(parser, message):
    parser.print_help()
    print("\n===== ERROR =======")
    print(message)
    print("See Usage Above!!!")
    sys.exit(False)
