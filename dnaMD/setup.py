#!/usr/bin/env python
#
#
# This file is part of do_x3dna
#
# Author: Rajendra Kumar
# Copyright (C) 2014  Rajendra Kumar
#
# do_x3dna uses 3DNA package (http://x3dna.org).
# Please cite the original publication of the 3DNA package:
# Xiang-Jun Lu & Wilma K. Olson (2003)
# 3DNA: a software package for the analysis, rebuilding and visualization of three-dimensional nucleic acid structures
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
#=============================================================================

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

def read(fname):
    return open(path.join(path.dirname(__file__), fname)).read()

setup(
    name = 'dnaMD',
    version = '2.0.6',

    # Required packages
    install_requires = [ 'numpy>=1.6',  'scipy>=0.9', 'h5py>=2.2.1' ],
    packages=find_packages(),

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={'console_scripts': [ 'dnaMD=dnaMD:main.main',], },

    # metadata for upload to pypi
    author = "Rajendra Kumar",
    author_email = "rjdkmr@gmail.com",
    url = 'https://github.com/rjdkmr/do_x3dna',
    description = "To analyze fluctuations in DNA or RNA structures in molecular dynamics trajectories",
    long_description = read('README.rst'),
    keywords = ["DNA", "MD Simulation", "Elastic Properties", "DNA Analysis MD Simulations"],
    license = 'GNU General Public License v3 (GPLv3)',
    classifiers = [
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: MacOS',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
