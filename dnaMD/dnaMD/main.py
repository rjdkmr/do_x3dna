
import os
import sys
import re


def main():

    options = {'vsTime':'Extract a parameter as a function of time.',
               'saveAsH5' : 'Store parameters from text file to HDF5 file.'
              }

    if len(sys.argv)<=1:
        show_help(options)
        sys.exit(-1)

    if sys.argv[1] not in options:
        print(' ERROR: "{0}" is not an accepted option.\n' .format(sys.argv[1]))
        show_help(options)
        sys.exit(-1)

    if sys.argv[1] == 'vsTime':
        from .commands import vsTime
        vsTime.main()

    if sys.argv[1] == 'saveAsH5':
        from .commands import saveAsH5
        saveAsH5.main()

def show_help(options):
    print(' ==============================')
    print(' Usage:')
    print(' dnaMD <Option>\n')
    print(' ---------------------')
    print(' Use following options:')
    print(' ---------------------\n')

    for tool in options:
        print('\t{0} : {1}\n'.format(tool, options[tool]))

    print(' ==============================')

if __name__ == '__main__':
    main()
