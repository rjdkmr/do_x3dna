
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


def main():

    options = {'vsTime':'Extract a parameter as a function of time.',
               'histogram': 'Parameter distribution during simulation.',
               'saveAsH5' : 'Store parameters from text file to HDF5 file.',
               'vsBPS' : 'Average parameters as a function of base-pair/step',
               'localDeformation' : 'Deformation of local parameters in probe DNA with respect to a reference DNA',
               'axisCurv' : 'Calculate global helical-axis, curvatures and tangents'
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

    if sys.argv[1] == 'histogram':
        from .commands import histogram
        histogram.main()

    if sys.argv[1] == 'vsBPS':
        from .commands import vsBPS
        vsBPS.main()

    if sys.argv[1] == 'localDeformation':
        from .commands import localDeformation
        localDeformation.main()

    if sys.argv[1] == 'saveAsH5':
        from .commands import saveAsH5
        saveAsH5.main()

    if sys.argv[1] == 'axisCurv':
        from .commands import axisCurv
        axisCurv.main()


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
