#!/usr/bin/env python
#
#
# This file is part of do_x3dna
#
# Author: Rajendra Kumar
# Copyright (C) 2014-2017  Rajendra Kumar
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

import numpy as np
import re
import os
import sys
import string
import random
import math
import subprocess as sub
try:
    from scipy.interpolate import splprep, splev
    from scipy.spatial import distance as sp_dist
    scipy_imported = True
except:
    scipy_imported = False

HAVE_H5PY = False
try:
    import h5py
    HAVE_H5PY = True
except:
    pass


## parameter count -- equal to basepairs
basePairParameters = ['shear', 'stretch', 'stagger', 'buckle', 'propeller', 'opening' ]
backboneDihedrals =  [  'alpha S1', 'beta S1', 'gamma S1', 'delta S1', 'epsilon S1', 'zeta S1', 'chi S1',
                        'alpha S2', 'beta S2', 'gamma S2', 'delta S2', 'epsilon S2', 'zeta S2', 'chi S2'
                     ]

## parameter count -- less than one of basepairs
baseStepParameters = [ 'shift', 'slide', 'rise', 'tilt', 'roll', 'twist' ]
helicalBaseStepParameters = [ 'x-disp', 'y-disp', 'h-rise', 'inclination', 'tip', 'h-twist' ]
helicalRadiusParameters = [ 'radius s-1', 'radius s-2']
helicalAxisParameters = [   'helical x-axis', 'helical y-axis', 'helical z-axis' ]
helicalAxisSmoothParameters =  [ 'helical x-axis smooth', 'helical y-axis smooth', 'helical z-axis smooth' ]
groovesParameters = [ 'major groove', 'major groove refined',
                      'minor groove', 'minor groove refined' ]
bendingParameters = [ 'curvature',  'tangent']



# Combine pasepair type parameters
basePairTypeParameters = basePairParameters + backboneDihedrals

# Combine basestep type parameters
baseStepTypeParameters = baseStepParameters + helicalBaseStepParameters \
                + helicalRadiusParameters + helicalAxisParameters \
                + helicalAxisSmoothParameters + groovesParameters \
                + bendingParameters


def getParameterType(param):
    if param in basePairTypeParameters:
        return 'bp'

    if param in baseStepTypeParameters:
        return 'bps'



class DNA:
    """DNA class stores all data obtained from the input files.

    **To initialize this class:** ::

        dna = DNA(60)       # 60 is the number of basepairs

    This class also contains several methods (functions) that are discussed in following sections.

    The data inside the class is organized as a dictionary in following architecture:

    .. code-block:: bash

                 |------- bp  ---- BP number      ----  Parameter array (1D/2D)
        data-----|
                 |------- bps ---- BP Step Number ----  Parameter array (1D/2D)


    Parameters
    ----------
    num_bp : int
        Number of base-pairs in the DNA.
    filename : str
        Name of HDF5 file
    startBP : int
        Number ID of first basepair.


    Attributes
    ----------
    num_bp : int
        Number of base-pairs in the DNA.
    num_step : int
        Number of base-steps in the DNA.
    filename : str
        Name of HDF5 file
    startBP : int
        Number ID of first basepair.
    smooth_axis : bool
        If axis is smoothened and global axis is already determined.
    data : dictionary
        Dictionary of data. All data are stored in this dictionary.
    time : ndarray
        Trajectory time.
    mask : ndarray
        Boolean array indicating which frame of trajectory should be masked.
    h5 : object
        h5py.File object to read the data directly from file.


    """

    def __init__(self, num_bp, filename=None, startBP=1):
        self.num_bp = num_bp
        self.num_step = num_bp - 1
        self.filename = filename
        self.startBP = startBP
        self.smooth_axis = False
        self.data = dict()
        self.data['bp'] = dict()
        self.data['bps'] = dict()
        self.time = []
        self.mask = None
        self.h5 = None

        # Check h5py is installed or not
        if self.filename is not None and not HAVE_H5PY:
            raise ImportError( 'h5py module is not installed. Please use "pip install h5py" or "pip3 install h5py" to install h5py.' )

        for i in range(self.startBP, self.startBP + num_bp):
            self.data['bp'][str(i)] = dict()

        for i in range(self.startBP, self.startBP + self.num_step):
            self.data['bps'][str(i)] = dict()


        self._openFile()

    def __del__(self):
        if self.filename is not None and self.h5 is not None:
            try:
                self.h5.close()
            except:
                pass

    def _openFile(self):
        """ Open the HDF5 file for reading/writing.

        This methad is called during the initialization of this class.
        """
        if self.filename is not None:
            self.h5 = h5py.File(self.filename)
        else:
            return

        if 'bp' not in self.h5:
            self.h5.create_group('bp')

        if 'bps' not in self.h5:
            self.h5.create_group('bps')

        for bp_type in ['bp', 'bps']:
            for bp_num in self.data[bp_type]:

                # Create group in case if new hdf5 file is opened
                if bp_num not in self.h5[bp_type]:
                    self.h5[bp_type].create_group(bp_num)

                # Sync data if old file is opened
                for parameters in self.h5[bp_type][bp_num]:
                    self.data[bp_type][bp_num][parameters] = self.h5[bp_type][bp_num][parameters]

        if 'mask' in self.h5:
            self.mask = self.h5['mask'][:]

        if 'time' in self.h5:
            self.time = self.h5['time'][:]

        if 'num_bp' in self.h5.attrs:
            self.num_bp = self.h5.attrs['num_bp']
        else:
            self.h5.attrs['num_bp'] = self.num_bp

        if 'num_step' in self.h5.attrs:
            self.num_step = self.h5.attrs['num_step']
        else:
            self.h5.attrs['num_step'] = self.num_step

        if 'startBP' in self.h5.attrs:
            self.startBP = self.h5.attrs['startBP']
        else:
            self.h5.attrs['startBP'] = self.startBP

    def _set_time(self, time):
        """ Set time in both class and hdf5 file
        """
        if len(self.time) == 0 :
            self.time = np.array(time)
            if self.h5 is not None:
                self.h5.create_dataset('time', self.time.shape, dtype=self.time.dtype, data=self.time, compression="gzip", shuffle=True, scaleoffset=3)
        else:
            if(len(time) != len(self.time)):
                raise AssertionError("\nTime or number of frame mismatch in input files.\n Exiting...\n")

    def _set_mask(self, mask):
        """ Set mask array in both class and hdf5 file
        """
        self.mask = mask.copy()
        if self.h5 is not None:
            if 'mask' in self.h5:
                self.h5.pop('mask')
            self.h5.create_dataset('mask', mask.shape, dtype=self.mask.dtype, data=mask, compression="gzip", shuffle=True)


    def _set_data(self, data, bp_type, bp_num, parameter, scaleoffset=3):
        if self.h5 is None:
            self.data[bp_type][bp_num][parameter] = data
        else:
            if parameter in self.h5[bp_type][bp_num]:          # Remove old data
                self.h5[bp_type][bp_num].pop( parameter )

            self.h5[bp_type][bp_num].create_dataset(parameter, data.shape, dtype=data.dtype, data=data, compression="lzf", shuffle=True, scaleoffset=scaleoffset)
            self.data[bp_type][bp_num][parameter] = self.h5[bp_type][bp_num][parameter]

    def get_parameters(self, parameter, bp, bp_range=True, masked=False):
        """To get the parameters over all frame for the given range of base pair/steps

        parameters
        ----------
        parameter : str
            Input parameter name.
            Parameter from following list are accepted:

                * ``shear``
                * ``stretch``
                * ``stagger``
                * ``buckle``
                * ``propeller``
                * ``opening``
                * ``shift``
                * ``slide``
                * ``rise``
                * ``tilt``
                * ``roll``
                * ``twist``
                * ``x-disp``
                * ``y-disp``
                * ``h-Rise``
                * ``inclination``
                * ``tip``
                * ``h-Twist``
                * ``helical X-axis``
                * ``helical Y-axis``
                * ``helical z-axis``
                * ``helical X-axis smooth``
                * ``helical Y-axis smooth``
                * ``helical z-axis smooth``
                * ``curvature``
                * ``tangent``
                * ``radius S-1``
                * ``radius S-2``
                * ``major Groove``
                * ``major Groove Refined``
                * ``minor Groove``
                * ``minor Groove Refined``
                * ``alpha S1``
                * ``beta S1``
                * ``gamma S1``
                * ``delta S1``
                * ``epsilon S1``
                * ``zeta S1``
                * ``chi S1``
                * ``alpha S2``
                * ``beta S2``
                * ``gamma S2``
                * ``delta S2``
                * ``epsilon S2``
                * ``zeta S2``
                * ``chi S2``

        bp : 1D list
            List of base-pairs to analyze
            Example: ::

                bp = [6]                                # bp_range = False
                bp = [4,15]                             # bp_range = True
                bp = range(4,15)                        # bp_range = False
                bp = np.arange(4,15)                    # bp_range = False
                bp = [2,5,6,7,9,12,18]                  # bp_range = False

        bp_range : bool
            ``Dfault=True``: As shown above, if ``True``, bp is taken as a range otherwise list or numpy array.

        masked : bool
            ``Dfault=False``: To skip specific frames/snapshots. dnaMD.DNA.mask array should
            be set to use this functionality. This array contains boolean (either ``True`` or ``False``)
            value for each frame to mask the frames. Presently, mask array is automatically generated
            during :meth:`DNA.generate_smooth_axis` method to skip those frames where 3D fitting
            curve was not successfull within the given critera.


        Returns
        -------
        parameters : 2D list
            ``parameters[bp][nframe] (2D list)``: where bp is number of base pairs/steps and nframe is total number of frames in the trajectory.

        """

        bpIndex, dum = get_idx_of_bp_parameters(bp, [], bp_range, startBP=self.startBP)
        append = False
        empty = False
        key = 'dummy'
        idx = 0
        data = []
        midx = []

        # Masking values according to mask array
        if masked and self.mask is None:
            print(" WARNING: mask is not set. Mask is set during helical axis smoothening. \n")
            masked = False

        for i in range(len(self.time)):
            if masked:
                if self.mask[i] == False:
                    midx.append(i)
            else:
                midx.append(i)

        param_type = getParameterType(parameter)
        if param_type is None:
            raise ValueError('ERROR: Incorrect parameter keyword: \"{0}\" .\n' .format(parameter))

        outData, bp_nums = [], []
        for i in range(len(bpIndex)):
            num = str( bpIndex[i]+self.startBP )
            if parameter in self.data[param_type][num]:
                tData = self.data[param_type][num][parameter][:]
                outData.append(tData[midx])
                bp_nums.append(num)
            else:
                raise ValueError('ERROR: The parameter \"{0}\" for base pair/step \"{1}\" is not set/loaded.\n' .format(parameter, num))

        return outData, bpIndex

    def time_vs_parameter(self, parameter, bp, merge=False, merge_method='mean', masked=False):
        """To get the parameter of either a specfic base-pair/step or a DNA segment as a function of time.

        parameters
        ----------
        parameter : str
            Name of a base-pair or base-step or helical parameter.
            For details about accepted keywords, see ``parameter`` in the
            method :meth:`DNA.get_parameters`.

        bp : 1D list or array
            base-pairs to analyze
            Example: ::

                bp = [6]                                  # merge = False
                bp = [4,15]                               # merge = True

        merge : bool
            ``Dfault=False``. As shown above, if ``True``, bp should a list of
            range otherwise a list of single value. If ``bp = True``, the
            parameter for the respective DNA segment could be merged or
            calculated by ``merge_method``.

        merge_method : str
            Method to calculate the parameter of a DNA segment from local
            parameters of all base-pairs/steps that are between the range
            given through ``bp``.
            Currently accepted keywords are as follows:

                * ``merge_method = mean``: Average of local parameters
                * ``merge_method = sum``: Sum of local parameters

        masked : bool
            ``Dfault=False``. To skip specific frames/snapshots.
            ``DNA.mask`` array should be set to use this functionality.
            This array contains boolean (either ``True`` or ``False``) value
            for each frame to mask the frames. Presently, mask array is
            automatically generated during :meth:`DNA.generate_smooth_axis` to
            skip those frames where 3D fitting curve was not successfull within
            the given critera.


        Returns
        -------
        time : 1D array
            Array containing time of each frame from trajectory
        value : 1D array
            Array containing parameter values for each frame from trajectory

        """
        if not (isinstance(bp, list) or isinstance(bp, np.ndarray)):
            raise AssertionError(
                "type %s is not list or np.ndarray" % type(bp))

        if (len(bp) > 1) and (merge == False):
            raise AssertionError(
                "bp %s contains more than two values, whereas merge=False. Use either one value in bp or merge=True" % bp)
            exit(1)

        if len(bp) == 1:
            merge = False

        if (merge == True) and not ((merge_method == 'mean') or (merge_method == 'sum')):
            raise AssertionError(
                "merge method %s is not available." % merge_method)

        # Masking values according to mask array
        midx = []
        if masked and self.mask is None:
            print(" WARNING: mask is not set. Mask is set during helical axis smoothening. \n")
            masked = False

        for i in range(len(self.time)):
            if masked:
                if self.mask[i] == False:
                    midx.append(i)
            else:
                midx.append(i)

        if len(bp) == 1:
            param_value, bp_idx = self.get_parameters(parameter, bp, bp_range=False, masked=masked)
        else:
            param_value, bp_idx = self.get_parameters(parameter, bp, bp_range=True, masked=masked)

        if (merge == True) and (merge_method == 'mean'):
            return self.time, np.mean(param_value, axis=0)

        elif (merge == True) and (merge_method == 'sum'):
            return self.time[midx], np.sum(param_value, axis=0)

        else:
            return self.time[midx], param_value[0]

    def parameter_distribution(self, parameter, bp, bins=30, merge=False, merge_method='mean', masked=False):
        """To get the parameter distribution of either a specfic base-pair/step or a DNA segment over the MD simulation.

        parameters
        ----------
        parameter : str
            Name of a base-pair or base-step or helical parameter
            For details about accepted keywords, see ``parameter`` in the method
            :meth:`DNA.get_parameters`.

        bp : 1D list or array
            base-pairs to analyze
            Example: ::

                bp = [6]                                  # merge = False
                bp = [4,15]                               # merge = True

        bins int
            Number of bins to calculate histogram

        merge : bool
            ``Dfault=False``: As shown above, if ``True``, bp should a list of
            range otherwise a list of single value. If ``bp = True``, the
            parameter for the respective DNA segment could be merged or
            calculated by ``merge_method``.

        merge_method : str
            Method to calculate the parameter of a DNA segment from local
            parameters of all base-pairs/steps that are between the range given
            through ``bp``.
            Currently accepted keywords are as follows:

                * ``merge_method = mean``: Average of local parameters
                * ``merge_method = sum``: Sum of local parameters

        masked : bool
            ``Dfault=False``. To skip specific frames/snapshots.
            ``DNA.mask`` array should be set to use this functionality.
            This array contains boolean (either ``True`` or ``False``) value
            for each frame to mask the frames. Presently, mask array is
            automatically generated during :meth:`DNA.generate_smooth_axis` to
            skip those frames where 3D fitting curve was not successfull within
            the given critera.


        Returns
        -------
        values : 1D array
            Array containing parameter values
        density : 1D array
            Array containing density for respective parameter values

        """
        if not (isinstance(bp, list) or isinstance(bp, np.ndarray)):
            raise AssertionError(
                "type %s is not list or np.ndarray" % type(bp))

        if (len(bp) > 1) and (merge == False):
            raise AssertionError(
                "bp %s contains more than two values, whereas merge=False. Use either one value in bp or merge=True" % bp)
            exit(1)

        if len(bp) == 1:
            merge = False

        if (merge == True) and not ((merge_method == 'mean') or (merge_method == 'sum')):
            raise AssertionError(
                "merge method %s is not available." % merge_method)
            exit(1)

        if len(bp) == 1:
            param_value, bp_idx = self.get_parameters(
                parameter, bp, bp_range=False, masked=masked)
        else:
            param_value, bp_idx = self.get_parameters(
                parameter, bp, bp_range=True, masked=masked)

        if (merge == True) and (merge_method == 'mean'):
            param_value = np.mean(param_value, axis=0)

        elif (merge == True) and (merge_method == 'sum'):
            param_value = np.sum(param_value, axis=0)

        else:
            param_value = param_value[0]

        density, bin_edges = np.histogram(param_value, bins=bins, density=True)
        bin_width = bin_edges[1] - bin_edges[0]

        density = np.insert(density, 0, 0.0)
        density = np.append(density, 0.0)

        values = []
        for i in range(len(bin_edges) - 1):
            values.append((bin_edges[i] + bin_edges[i + 1]) / 2)

        values = np.asarray(values)
        values = np.append(values, values[-1] + bin_width)
        values = np.insert(values, 0, values[0] - bin_width)

        return np.array(values), density

    def set_base_pair_parameters(self, filename, bp, parameters=['shear', 'stretch', 'stagger', 'buckle', 'propeller', 'opening'], bp_range=True):
        """To read and store basepairs parameters (shear, stretch, stagger, buckle, propeller and opening) from an input file.

        parameters
        ----------
        filename : str
            Input file, which is generated from do_x3dna. e.g. L-BP_g.dat

        bp : 1D list or array
            base-pairs to analyze
            Example: ::

                bp = [6]                                # bp_range = False
                bp = [4,15]                             # bp_range = True
                bp = range(4,15)                        # bp_range = False
                bp = np.arange(4,15)                    # bp_range = False
                bp = [2,5,6,7,9,12,18]                  # bp_range = False

        parameters : 1D list
            List of base-pairs parameters as follows:

                * ``shear``
                * ``stretch``
                * ``stagger``
                * ``buckle4``
                * ``propeller``
                * ``opening6``

            By default all six parameters will be extracted from the file.

        bp_range : bool
            ``Dfault=True``: As shown above, if ``True``, bp is taken as a range otherwise list or numpy array.

        """
        if not (isinstance(bp, list) or isinstance(bp, np.ndarray)):
            raise AssertionError(
                "type %s is not list or np.ndarray" % type(bp))
        if not (isinstance(parameters, list) or isinstance(parameters, np.ndarray)):
            raise AssertionError(
                "type %s is not list or np.ndarray" % type(parameters))

        targetParameters = { 1:'shear', 2:'stretch', 3:'stagger', 4:'buckle', 5:'propeller', 6:'opening' }
        targetParametersReverse = dict((v,k) for k,v in targetParameters.items())

        # Check if requested parameters found within input file
        gotParametersInputFile = checkParametersInputFile(filename)
        if gotParametersInputFile is None:
            raise IOError(' Something wrong in input file {0}.\n Cannot read parameters.\n File should be an ouput from do_x3dna.'.format(filename))
        for parameter in parameters:
            if parameter not in gotParametersInputFile:
                raise ValueError(' Parameter {0} not found in input file. \n This file contains following parameters: \n {1}'.format(parameter, gotParametersInputFile))


        InputParamIndex = []
        for parameter in parameters:
            if parameter in targetParameters.values():
                InputParamIndex.append( targetParametersReverse[parameter] )
            else:
                print('\nWARNING: base pair parameters \"{0}\" not accepted. Skipping it !!\n' .format(parameter))

        if not InputParamIndex:
            raise ValueError("No acceptable base-pair parameters found!!!")

        data, time = read_param_file(filename, InputParamIndex, bp, bp_range, startBP=self.startBP)
        self._set_time(time)
        bpIndex, OutParamIndex = get_idx_of_bp_parameters(bp, InputParamIndex, bp_range, startBP=self.startBP)

        for i in range(len(data)):
            for j in range(len(data[i])):
                bp_num = str( bpIndex[i]+self.startBP )
                param = targetParameters[OutParamIndex[j]+1]
                self._set_data(data[i][j], 'bp', bp_num, param, scaleoffset=2)

    def set_major_minor_groove(self, filename, bp_step, parameters=[  'minor groove', 'minor groove refined', 'major groove', 'major groove refined'], step_range=True):
        """To read and store Major and Minor grooves from an input file.

        * Minor groove : direct P-P distance
        * Minor Grrove Refined : refined P-P distance which take into account the directions of the sugar-phosphate backbones
        * Major groove : direct P-P distance
        * Major Grrove Refined : refined P-P distance which take into account the directions of the sugar-phosphate backbones

        .. warning::

                * The major and minor grooves (direct P-P) cannot be calculated for first and last two base-steps
                * The major and minor grooves (refined P-P) cannot be calculated for first and last three base-steps

        Parameters
        ----------
        filename : str
            Input file, which is generated from do_x3dna. e.g. MGroove_g.dat

        bp_step : 1D list or array
            base-steps to analyze
            Example: ::

                bp_step = [6]                                # step_range = False
                bp_step = [4,15]                             # step_range = True
                bp_step = range(4,15)                        # step_range = False
                bp_step = np.arange(4,15)                    # step_range = False
                bp_step = [2,5,6,7,9,12,18]                  # step_range = False

        parameters : 1D list
            List of groove parameter names as follows:
                * ``minor groove``
                * ``minor grrove refined``
                * ``major groove``
                * ``major grrove refined``

            By default all four groove parameters will be extracted from the file.

        step_range : bool
            ``Dfault=True``: As shown above, if ``True``, bp_step is taken as a range otherwise list or numpy array.

        """
        if not (isinstance(bp_step, list) or isinstance(bp_step, np.ndarray)):
            raise AssertionError(
                "type %s is not list or np.ndarray" % type(bp_step))

        if not (isinstance(parameters, list) or isinstance(parameters, np.ndarray)):
            raise AssertionError(
                "type %s is not list or np.ndarray" % type(parameters))

        targetParameters = { 1:'minor groove', 2:'minor groove refined', 3:'major groove', 4:'major groove refined' }
        targetParametersReverse = dict((v,k) for k,v in targetParameters.items())

        # Check if requested parameters found within input file
        gotParametersInputFile = checkParametersInputFile(filename)
        if gotParametersInputFile is None:
            raise IOError(' Something wrong in input file {0}.\n Cannot read parameters.\n File should be an ouput from do_x3dna.'.format(filename))
        for parameter in parameters:
            if parameter not in gotParametersInputFile:
                raise ValueError(' Parameter {0} not found in input file. \n This file contains following parameters: \n {1}'.format(parameter, gotParametersInputFile))

        InputParamIndex = []
        for parameter in parameters:
            if parameter in targetParameters.values():
                InputParamIndex.append( targetParametersReverse[parameter] )
            else:
                print('\nWARNING: base pair parameters \"{0}\" not accepted. Skiiping it !!\n' .format(parameter))

        data, time = read_param_file(filename, InputParamIndex, bp_step, step_range, word=True, startBP=self.startBP)
        self._set_time(time)

        bpIndex, OutParamIndex = get_idx_of_bp_parameters(bp_step, InputParamIndex, step_range, startBP=self.startBP)

        for i in range(len(data)):
            for j in range(len(data[i])):
                # terminal base-steps do not have major and minor grooves
                if data[i][j][0] != None:
                    bp_num = str( bpIndex[i]+self.startBP )
                    param = targetParameters[OutParamIndex[j]+1]
                    self._set_data(np.asarray(data[i][j], dtype=np.float), 'bps', bp_num, param, scaleoffset=1)

    def set_backbone_dihedrals(self, filename, bp, parameters='All', bp_range=True):
        """To read and store backbone dihedrals (alpha, beta, gamma, delta, epsilon and zeta) and chi dihedral of both strands from an input file.

        .. note::

                * alpha:   O3'(i-1)-P-O5'-C5'
                * beta:    P-O5'-C5'-C4'
                * gamma:   O5'-C5'-C4'-C3'
                * delta:   C5'-C4'-C3'-O3'
                * epsilon: C4'-C3'-O3'-P(i+1)
                * zeta:    C3'-O3'-P(i+1)-O5'(i+1)
                * chi for pyrimidines(Y): O4'-C1'-N1-C2
                * chi for purines(R): O4'-C1'-N9-C4


        Parameters
        ----------
        filename : str
            Input file, which is generated from do_x3dna. e.g. BackBoneCHiDihedrals_g.dat

        bp : 1D list or array
            base-pairs to analyze
            Example: ::

                bp = [6]                                # bp_range = False
                bp = [4,15]                             # bp_range = True
                bp = range(4,15)                        # bp_range = False
                bp = np.arange(4,15)                    # bp_range = False
                bp = [2,5,6,7,9,12,18]                  # bp_range = False

        parameters : str or 1D list
            Either ``All`` to extract all angles or list of angles as follows:

                * ``alpha S1``
                * ``beta S1``
                * ``gamma S1``
                * ``delta S1``
                * ``epsilon S1``
                * ``zeta S1``
                * ``chi S1``
                * ``alpha S2``
                * ``beta S2``
                * ``gamma S2``
                * ``delta S2``
                * ``epsilon S2``
                * ``zeta S2``
                * ``chi S2``

        bp_range : bool
            ``Dfault=True``: As shown above, if ``True``, bp is taken as a range otherwise list or numpy array.

        """
        if not (isinstance(bp, list) or isinstance(bp, np.ndarray)):
            raise AssertionError(
                "type %s is not list or np.ndarray" % type(bp))

        if not (isinstance(parameters, list) or isinstance(parameters, np.ndarray)):
            if parameters == 'All':
                parameters = [ 'alpha S1', 'beta S1', 'gamma S1', 'delta S1', 'epsilon S1', 'zeta S1', 'chi S1',
                               'alpha S2', 'beta S2', 'gamma S2', 'delta S2', 'epsilon S2', 'zeta S2', 'chi S2'  ]
            else:
                raise ValueError(" ERROR: {0} is not accepted parameters!!! It should be either \"All\" or list of parameter names.".format(parameters) )

        targetParameters = { 1:'alpha S1', 2:'beta S1',  3:'gamma S1',  4:'delta S1',  5:'epsilon S1',  6:'zeta S1',  7:'chi S1',
                             8:'alpha S2', 9:'beta S2', 10:'gamma S2', 11:'delta S2', 12:'epsilon S2', 13:'zeta S2', 14:'chi S2'  }
        targetParametersReverse = dict((v,k) for k,v in targetParameters.items())

        # Check if requested parameters found within input file
        gotParametersInputFile = checkParametersInputFile(filename)
        if gotParametersInputFile is None:
            raise IOError(' Something wrong in input file {0}.\n Cannot read parameters.\n File should be an ouput from do_x3dna.'.format(filename))
        for parameter in parameters:
            if parameter not in gotParametersInputFile:
                raise ValueError(' Parameter {0} not found in input file. \n This file contains following parameters: \n {1}'.format(parameter, gotParametersInputFile))

        InputParamIndex = []
        for parameter in parameters:
            if parameter in targetParameters.values():
                InputParamIndex.append( targetParametersReverse[parameter] )
            else:
                print('\nWARNING: base pair parameters \"{0}\" not accepted. Skipping it !!\n' .format(parameter))

        if not InputParamIndex:
            raise ValueError("No acceptable base-pair parameters found!!!")

        data, time = read_param_file(filename, InputParamIndex, bp, bp_range, word=True, startBP=self.startBP)
        self._set_time(time)

        bpIndex, OutParamIndex = get_idx_of_bp_parameters(bp, InputParamIndex, bp_range, startBP=self.startBP)

        for i in range(len(data)):
            bp_num = str( bpIndex[i]+self.startBP )
            for j in range(len(data[i])):
                if data[i][j][0] != None:
                    param = targetParameters[OutParamIndex[j]+1]
                    self._set_data(np.asarray(data[i][j], dtype=np.float), 'bp', bp_num, param, scaleoffset=1)

    def set_helical_radius(self, filename, bp, atomname='P', full=False, bp_range=True):
        """To read and set local helical radius of both strand

        Parameters
        ----------
        filename : str
            Input file, which is generated from do_x3dna. e.g. HelixRad_g.dat

        bp : 1D list or array
            base-pairs to analyze
            Example: ::

                bp = [6]                                # bp_range = False
                bp = [4,15]                             # bp_range = True
                bp = range(4,15)                        # bp_range = False
                bp = np.arange(4,15)                    # bp_range = False
                bp = [2,5,6,7,9,12,18]                  # bp_range = False

        atomname : str
            Atom name to consider for the DNA helix (accepted keywords:
            ``P``, ``O4*``, ``O4'``, ``C1*`` and ``C1``)

        full : bool
            To calculate full helical radius. Overrides atomname option and
            uses atom ``P``, and 1 A is added to the radius calculated by 3DNA
            package

        bp_range : bool
            ``Dfault=True``: As shown above, if ``True``, bp is taken as a
            range otherwise list or numpy array.

        """
        if not (isinstance(bp, list) or isinstance(bp, np.ndarray)):
            raise AssertionError(
                "type %s is not list or np.ndarray" % type(bp))

        if not ((atomname == 'P') or (atomname == 'O4*') or (atomname == 'C1*') or (atomname == 'O4\'') or (atomname == 'C1\'')):
            print (
                '\n This atomname {0} is not implemented... Exiting\n' .format(atomname))

        # Check if requested parameters found within input file
        gotParametersInputFile = checkParametersInputFile(filename)
        if gotParametersInputFile is None:
            raise IOError(' Something wrong in input file {0}.\n Cannot read parameters.\n File should be an ouput from do_x3dna.'.format(filename))
        for p in helicalRadiusParameters:
            if p not in gotParametersInputFile:
                raise ValueError(' Parameter {0} not found in input file. \n This file contains following parameters: \n {1}'.format(parameter, gotParametersInputFile))

        parameter = []
        if (atomname == 'P') or (full):
            parameter = [1, 4]

        if (atomname == 'O4*' or atomname == 'O4\'') and (not full):
            parameter = [2, 5]

        if (atomname == 'C1*' or atomname == 'C1\'') and (not full):
            parameter = [3, 6]


        data, time = read_param_file(filename, [1, 2, 3, 4, 5, 6], bp, bp_range, startBP=self.startBP)
        self._set_time(time)

        bp_idx, param_idx = get_idx_of_bp_parameters(bp, parameter, bp_range, startBP=self.startBP)

        if full:
            data = np.add(data, 1.0)

        for i in range(len(data)):
            bp_num = str( bp_idx[i]+self.startBP )
            if (atomname == 'P') or (full):
                self._set_data(data[i][0], 'bps', bp_num, 'radius s-1', scaleoffset=1)
                self._set_data(data[i][3], 'bps', bp_num, 'radius s-2', scaleoffset=1)

            if (atomname == 'O4*' or atomname == 'O4\''):
                self._set_data(data[i][1], 'bps', bp_num, 'radius s-1', scaleoffset=1)
                self._set_data(data[i][4], 'bps', bp_num, 'radius s-2', scaleoffset=1)

            if (atomname == 'C1*' or atomname == 'C1\''):
                self._set_data(data[i][2], 'bps', bp_num, 'radius s-1', scaleoffset=1)
                self._set_data(data[i][5], 'bps', bp_num, 'radius s-2', scaleoffset=1)

    def set_base_step_parameters(self, filename, bp_step, parameters='All', step_range=True, helical=False):
        """To read and store base-step (Shift, Slide, Rise, Tilt, Roll and Twist) and helical base-step (X-disp, Y-disp, h-Rise, Inclination, Tip and h-Twist) parameters from an input file

        Parameters
        ----------
        filename : str
            Input file, which is generated from do_x3dna. e.g. L-BPS_g.dat or L-BPH_g.dat.

        bp_step : 1D list or array
            base-steps to analyze
            Example: ::

                bp_step = [6]                                # step_range = False
                bp_step = [4,15]                             # step_range = True
                bp_step = range(4,15)                        # step_range = False
                bp_step = np.arange(4,15)                    # step_range = False
                bp_step = [2,5,6,7,9,12,18]                  # step_range = False

        parameters : str or 1D list
            Either ``All`` to extract all parameter available in input file
            or list of either base-step or helical base-step parameter as follows:

            If helical = ``False``:

                * ``shift``
                * ``slide``
                * ``rise``
                * ``tilt``
                * ``roll``
                * ``twist``

            If helical = ``True``:

                * ``X-disp``
                * ``Y-disp``
                * ``h-Rise``
                * ``Inclination``
                * ``Tip``
                * ``h-Twist``

        step_range : bool
            ``Dfault=True``: As shown above, if ``True``, bp_step is taken as a range otherwise list or numpy array.

        helical : bool
            If ``True``, parameters in input file will be considered as helical base-steps parameters
            If ``False``, parameters will be considered as base-steps parameters.

        """
        if not (isinstance(bp_step, list) or isinstance(bp_step, np.ndarray)):
            raise AssertionError(
                "type %s is not list or np.ndarray" % type(bp_step))

        if not (isinstance(parameters, list) or isinstance(parameters, np.ndarray)):
            if parameters == 'All' and not helical:
                parameters = ['shift', 'slide', 'rise', 'tilt', 'roll', 'twist']
            elif parameters == 'All' and helical:
                parameters = ['x-disp', 'y-disp', 'h-rise', 'inclination', 'tip', 'h-twist']
            else:
                raise ValueError(" ERROR: {0} is not accepted parameters!!! It should be either \"All\" or list of parameter names.".format(parameters) )

        targetParameters = None
        targetParametersReverse = None
        if helical:
            targetParameters = { 1:'x-disp', 2:'y-disp', 3:'h-rise', 4:'inclination', 5:'tip', 6:'h-twist' }
        else:
            targetParameters = { 1:'shift', 2:'slide', 3:'rise', 4:'tilt', 5:'roll', 6:'twist' }
        targetParametersReverse = dict((v,k) for k,v in targetParameters.items())

        # Check if requested parameters found within input file
        gotParametersInputFile = checkParametersInputFile(filename)
        if gotParametersInputFile is None:
            raise IOError(' Something wrong in input file {0}.\n Cannot read parameters.\n File should be an ouput from do_x3dna.'.format(filename))
        for parameter in parameters:
            if parameter not in gotParametersInputFile:
                raise ValueError(' Parameter {0} not found in input file. \n This file contains following parameters: \n {1}'.format(parameter, gotParametersInputFile))

        InputParamIndex = []
        for parameter in parameters:
            if parameter in targetParameters.values():
                InputParamIndex.append( targetParametersReverse[parameter] )
            else:
                print('\nWARNING: base pair parameters \"{0}\" not accepted. Skiiping it !!\n' .format(parameter))

        if not InputParamIndex:
            raise ValueError("No acceptable base-pair parameters found!!!")

        data, time = read_param_file(filename, InputParamIndex, bp_step, step_range, startBP=self.startBP)
        self._set_time(time)

        bpIndex, OutParamIndex = get_idx_of_bp_parameters(bp_step, InputParamIndex, step_range, startBP=self.startBP)

        for i in range(len(data)):
            for j in range(len(data[i])):
                bp_num = str( bpIndex[i]+self.startBP )
                param = targetParameters[OutParamIndex[j]+1]
                self._set_data(data[i][j], 'bps', bp_num, param, scaleoffset=2)

    def get_mean_error(self, bp, parameter, err_type='std', bp_range=True, merge_bp=1, merge_method='mean', masked=False, tool='g_analyze'):
        """To calculate average and error of the given parameter for the gieven set of base-pairs/steps

        .. warning::
                To calculate errors by using ``error = 'acf'`` or ``error = 'block'``,
                GROMACS tool ``g_analyze`` or ``gmx analyze`` should be present in ``$PATH``.

        Parameters
        ----------

        bp : 1D list or array
            base-pairs to analyze
            Example: ::

                bp = [6]                                # bp_range = False
                bp = [4,15]                             # bp_range = True
                bp = range(4,15)                        # bp_range = False
                bp = np.arange(4,15)                    # bp_range = False
                bp = [2,5,6,7,9,12,18]                  # bp_range = False

        parameter : str
            Name of a parameter. For details about accepted keywords, see ``parameter`` in the method :meth:`DNA.get_parameters`.

        error : str
            Method of error estimation.
            Currently accepted method as follows:

                * ``error = 'std'``   : Standard Deviation
                * ``error = 'acf'``   : Standard error using autocorrelation time (requires: ``g_analyze`` or ``gmx analyze``)
                * ``error = 'block'`` : Standard error using block averaging method (requires: ``g_analyze`` or ``gmx analyze``)

        bp_range : bool
            ``Dfault=True``: As shown above, if ``True``, bp is taken as a range otherwise list or numpy array.

        merge_bp : int
            Number of base-pairs or steps to merge for creating the small DNA segments

        merge_method : str
            Method to calculate the parameter of a DNA segment from local
            parameters of all base-pairs/steps that are between the range
            given through ``bp``.
            Currently accepted keywords are as follows:

                * ``merge_method = mean``: Average of local parameters
                * ``merge_method = sum``: Sum of local parameters

        masked : bool
            ``Dfault=False``. To skip specific frames/snapshots.
            ``DNA.mask`` array should be set to use this functionality.
            This array contains boolean (either ``True`` or ``False``) value
            for each frame to mask the frames. Presently, mask array is
            automatically generated during :meth:`DNA.generate_smooth_axis` to
            skip those frames where 3D fitting curve was not successfull within
            the given critera.

        tool : str
            Gromacs tool ``g_analyze`` or ``gmx analyze`` or ``gmx_mpi analyze`` etc.
            It will be used to calculate autocorrelation time or  block averaging error.
            It should be present in ``$PATH``

        Returns
        -------
        basepairs/steps : 1D array
            Number of base pair-steps. If ``merge_bp>1``, middle number will be returned.
        values : 1D array
            Average values of the parameter
        errors : 1D array
            Error values for corresponding average values

        """

        merge = False

        if(not bp_range) and (merge_bp > 1):
            print (
                "\nERROR: Merging of base pairs/steps only possible with given base pair/steps range\n")
            exit(1)

        if(bp_range) and (merge_bp > 1):
            merge = True

        if (merge == True) and not ((merge_method == 'mean') or (merge_method == 'sum')):
            raise AssertionError(
                "merge method %s is not available." % merge_method)
            exit(1)

        data, bp_idx = self.get_parameters(parameter, bp, bp_range, masked=masked)

        bp_number = np.add(bp_idx, self.startBP)
        data = np.array(data)

        # Merging base pairs/step data for the given parameters
        merge_data = []
        mid_bin = []
        if(merge):
            i = 0
            while(1):
                start = i
                if((i + merge_bp) >= len(bp_idx)):
                    end = i + (len(bp_idx) - i)
                else:
                    end = i + merge_bp

                if(start < end):
                    mid_bin.append(
                        int(start + (end - start) / 2) + bp_number[0])
                    if (merge_method == 'mean'):
                        merge_data.append(np.mean(data[start:end], axis=0))
                    if (merge_method == 'sum'):
                        merge_data.append(np.sum(data[start:end], axis=0))

                if(i >= len(bp_idx)):
                    break

                i += merge_bp

            merge_data = np.array(merge_data)

        if (err_type == 'std'):
            if(merge):
                error = np.std(merge_data, axis=1)
            else:
                error = np.std(data, axis=1)

        if (err_type == 'acf') or (err_type == 'block'):
            if(merge):
                error = get_error(self.time, merge_data, len(mid_bin), err_type=err_type, tool=tool)
            else:
                error = get_error(self.time, data, len(bp_idx), err_type=err_type, tool=tool)

        if(merge):
            return mid_bin, np.mean(merge_data, axis=1), error
        else:
            return bp_number, np.mean(data, axis=1), error

    def set_helical_axis(self, filename, step_range=False, step=None):
        """
        To read and set local helical-axis postions from an input file.

        Parameters
        ----------
        filename : str
            Input file, which is generated from do_x3dna. e.g. HelAxis_g.dat
        step_range : bool
            * ``step_range = True`` : read axis coordinates of base-steps for the given range of base-steps
            * ``step_range = False``: read axis coordinates of all base-steps
        step : list
            List containing lower and higher limit of base-steps range.
                * This option only works with ``step_range=True``.
                * This list should not contain more than two number.
                * First number should be less than second number.

            Example for base-step 4 to 15:
                ``step = [4,15]         # step_range = True``

        """

        if (step_range):
            if not isinstance(step, list):
                raise AssertionError("type %s is not list" % type(step))
            if (len(step) > 2):
                print (
                    "ERROR: Range for helical axis should be list of two numbers, e.g. step=[1, 20] \n")
                exit(1)

        if (step_range) and (step == None):
            raise ValueError(
                "See, documentation for step  and step_range usage!!!")

        # Check if requested parameters found within input file
        gotParametersInputFile = checkParametersInputFile(filename)
        if gotParametersInputFile is None:
            raise IOError(' Something wrong in input file {0}.\n Cannot read parameters.\n File should be an ouput from do_x3dna.'.format(filename))
        for p in helicalAxisParameters:
            if p not in gotParametersInputFile:
                raise ValueError(' Parameter {0} not found in input file. \n This file contains following parameters: \n {1}'.format(parameter, gotParametersInputFile))

        targetParameters = { 1:'helical x-axis', 2:'helical y-axis', 3:'helical z-axis' }

        if (step_range):
            if (len(step) != 2):
                raise ValueError("See, documentation for step usage!!!")

            if step[0] > step[1]:
                raise ValueError("See, documentation for step usage!!!")
            data, time = read_param_file(filename, [1, 2, 3], step, True, startBP=self.startBP)
        else:
            data, time = read_param_file(filename, [1, 2, 3], [1, self.num_step], True, startBP=self.startBP)

        self._set_time(time)

        if (step_range):
            bp_idx, param_idx = get_idx_of_bp_parameters(step, [], True, startBP=self.startBP)
        else:
            bp_idx, param_idx = get_idx_of_bp_parameters([1, self.num_step], [], True, startBP=self.startBP)

        for i in range(len(data)):
            for j in range(len(data[i])):
                bp_num = str( bp_idx[i]+self.startBP )
                param = targetParameters[j+1]
                self._set_data(data[i][j], 'bps', bp_num, param, scaleoffset=2)

    def generate_smooth_axis(self, step_range=False, step=None, smooth=500.0, spline=3, fill_point=6, cut_off_angle=20):
        """To determine the global helical axis by smoothening local axis using spline interpolation.

        .. note::
            A 3D curve is fitted on local helical axis that are calculated using
            ``do_x3dna`` tool. Sometimes in few frames, fitting **may not** be
            accurate and produces artifact. To record these frames, ``DNA.mask``
            array containing boolean values are generated. If value is ``True``,
            fitting might not be correct and vice versa. This array could be
            used in later analysis to skip/mask the frames containing inaccurate
            axis.

        .. warning::
                This function requires `SciPy package <http://www.scipy.org/>`_.

        Parameters
        ----------
        step_range : bool
            * ``step_range = True`` : Smoothen axis for the given range of base-steps
            * ``step_range = False``: Smoothen axis for entire DNA. If original helical-axis of any base-step will be found to be not available, error will be raised.

        step : list
            List containing lower and higher limit of base-steps range.
                * This option only works with ``step_range=True``.
                * This list should not contain more than two number.
                * First number should be less than second number.

            Example for base-step 4 to 15:
                ``step = [4,15]         # step_range = True``

        smooth : float
            A smoothing condition. For more details, see about ``s = None``, which is paased into
            `scipy.interpolate.splprep() <https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.splprep.html>`_ method.

            .. warning ::
                * Lower value may lead to an artifact of local sharp kink in the smoothed axis.
                * Higher value may lead to the calculation of wrong helical axis.

        spline : int
            Degree of spline. For more details, see about ``k = 3``, which is paased into
            `scipy.interpolate.splprep() <https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.splprep.html>`_ method.

        fill_point : int
            Number of intrapolated points between two adjacent helical-axis coordinates.

        cut_off_angle : float
            Cut-off bending angle to define sharp kink in fitted curve. If angle
            in fitted curve is larger than this cut-off, refitting will be performed
            after deleting few of the original helical axis positions. If after
            this deletions, bending angle will not reduce below cut-off angle,
            value of ``smooth`` will be increased by 100 and entire cycle of
            fitting-refitting will be performed. When, value of ``smooth`` increases
            to more than 10000 during this fitting-refitting cycles, fitting process
            will be stopped with a warning message.

            To record the frames with bad fitting, ``True`` value will be
            stored in ``DNA.mask`` array for respective frame.

        """

        if not scipy_imported:
            raise ImportError(
                "SciPy package is not available. Please visit http://www.scipy.org/install.html for download and installation instructions.\n")

        if (step_range) and (step == None):
            raise ValueError(
                "See, documentation for step  and step_range usage!!!")

        bp_idx = []

        if step_range:
            if (len(step) != 2):
                raise ValueError("See, documentation for step usage!!!")

            if step[0] > step[1]:
                raise ValueError("See, documentation for step usage!!!")
            RawX, bp_idx = self.get_parameters( 'helical x-axis', step, bp_range=True)
            RawY, dummy = self.get_parameters( 'helical y-axis', step, bp_range=True)
            RawZ, dummy = self.get_parameters( 'helical z-axis', step, bp_range=True)
        else:
            RawX, bp_idx = self.get_parameters( 'helical x-axis', [1, self.num_step], bp_range=True)
            RawY, dummy = self.get_parameters( 'helical y-axis', [1, self.num_step], bp_range=True)
            RawZ, dummy = self.get_parameters( 'helical z-axis', [1, self.num_step], bp_range=True)

        RawX = np.array(RawX).T
        RawY = np.array(RawY).T
        RawZ = np.array(RawZ).T

        smoothX, smoothY, smoothZ = [], [], []

        nframes = len(self.time)

        maskArray = np.zeros(nframes, dtype=bool)

        for i in range(nframes):

            frame_number = i + 1
            if((frame_number < 100) and (frame_number % 10 == 0)) or ((frame_number < 1000) and (frame_number % 100 == 0)) or (frame_number % 1000 == 0):
                sys.stdout.write("\rFitting spline curve on helcial axis of frame %d out of %d frames" % (
                    frame_number, nframes))
                sys.stdout.flush()

            xsmooth, ysmooth, zsmooth, mask = fit_axis(bp_idx, frame_number, RawX[i], RawY[i], RawZ[i], smooth, spline, fill_point, cut_off_angle)
            maskArray[i] = mask

            smoothX.append(xsmooth)
            smoothY.append(ysmooth)
            smoothZ.append(zsmooth)

        sys.stdout.write("\nFinished spline curve fitting...\n")
        sys.stdout.flush()

        smoothX = np.asarray(smoothX).T
        smoothY = np.asarray(smoothY).T
        smoothZ = np.asarray(smoothZ).T

        self.smooth_axis = True
        for i in range(len(bp_idx)):
            bp_num = str( bp_idx[i]+self.startBP )
            self._set_data(smoothX[i], 'bps', bp_num, 'helical x-axis smooth', scaleoffset=2)
            self._set_data(smoothY[i], 'bps', bp_num, 'helical y-axis smooth', scaleoffset=2)
            self._set_data(smoothZ[i], 'bps', bp_num, 'helical z-axis smooth', scaleoffset=2)

        # Set mask array
        self._set_mask(maskArray)

    def write_haxis_pdb(self, filename='helical_axis.pdb', step_range=False, step=None, write_smooth_axis=True, write_orig_axis=False, write_curv=False, scale_curv=1):
        """To write trajectory of helcial-axis as a PDB format file.

        Both local helical axis and global (smoothened) axis can be written to PDB file.
        For global axis, curvature could be written in B-factor field of PDB file.

        Parameters
        ----------
        filename : str
            Name of the output PDB format file.

        step_range : bool
            * ``step_range = True`` : Smoothen axis for the given range of base-steps
            * ``step_range = False``: Smoothen axis for entire DNA. If original helical-axis of any base-step will be found to be not available, error will be raised.

        step : list
            List containing lower and higher limit of base-steps range.
                * This option only works with ``step_range=True``.
                * This list should not contain more than two number.
                * First number should be less than second number.

            Example for base-step 4 to 15:
                ``step = [4,15]         # step_range = True``

        write_smooth_axis : bool
            Write coordinates of smoothed helical axis as chain A.

        write_orig_axis : bool
            Write coordinates of original helical axis (output from do_x3dna) as chain B.

        write_curv : bool
            Write curvature of smoothed helical axis in B-factor coloumn of PDB file.

        scale_curv : int
            Scaling of curvature. ``curvature * scale_curv`` is written in  B-factor coloumn of PDB file.

        """

        if (step_range) and (step == None):
            raise ValueError(
                "See, documentation for step  and step_range usage!!!")

        if not write_orig_axis and not write_smooth_axis:
            raise ValueError(
                "Nothing to write as both \"write_orig_axis=Flase\" and \"write_smooth_axis=False\" !!!")

        if step_range:
            if (len(step) != 2):
                raise ValueError("See, documentation for step usage!!!")

            if step[0] > step[1]:
                raise ValueError("See, documentation for step usage!!!")

            # Orignal helical axis
            if (write_orig_axis):
                RawX, bp_idx = self.get_parameters(
                    'helical x-axis', step, bp_range=True)
                RawY, dummy = self.get_parameters(
                    'helical y-axis', step, bp_range=True)
                RawZ, dummy = self.get_parameters(
                    'helical z-axis', step, bp_range=True)

            # Smoothed helical axis
            if (write_smooth_axis):
                SmoothX, bp_idx = self.get_parameters(
                    'helical x-axis smooth', step, bp_range=True)
                SmoothY, bp_idx = self.get_parameters(
                    'helical y-axis smooth', step, bp_range=True)
                SmoothZ, bp_idx = self.get_parameters(
                    'helical z-axis smooth', step, bp_range=True)

            # curvature
            if (write_curv):
                curvature, bp_idx = self.get_parameters(
                    'curvature', step, bp_range=True)

        else:

            # Orignal helical axis
            if (write_orig_axis):
                RawX, bp_idx = self.get_parameters(
                    'helical x-axis', [1, self.num_step], bp_range=True)
                RawY, dummy = self.get_parameters(
                    'helical y-axis', [1, self.num_step], bp_range=True)
                RawZ, dummy = self.get_parameters(
                    'helical z-axis', [1, self.num_step], bp_range=True)

            # Smoothed helical axis
            if (write_smooth_axis):
                SmoothX, bp_idx = self.get_parameters(
                    'helical x-axis smooth', [1, self.num_step], bp_range=True)
                SmoothY, bp_idx = self.get_parameters(
                    'helical y-axis smooth', [1, self.num_step], bp_range=True)
                SmoothZ, bp_idx = self.get_parameters(
                    'helical z-axis smooth', [1, self.num_step], bp_range=True)

            # curvature
            if (write_curv):
                curvature, bp_idx = self.get_parameters(
                    'curvature', [1, self.num_step], bp_range=True)

        if (write_orig_axis):
            RawX = np.array(RawX).T
            RawY = np.array(RawY).T
            RawZ = np.array(RawZ).T

        if (write_smooth_axis):
            SmoothX = np.array(SmoothX).T
            SmoothY = np.array(SmoothY).T
            SmoothZ = np.array(SmoothZ).T

        if (write_curv):
            curvature = np.array(curvature).T

        f = open(filename, 'w')

        for i in range(len(self.time)):
            f.write('%-6s    %4d\n' % ("MODEL", i + 1))

            bfactor = 0.00

            if (write_smooth_axis):
                for j in range(len(SmoothX[i])):

                    if (write_curv):
                        bfactor = curvature[i][j] * scale_curv

                    f.write('%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n' % ("ATOM", j + 1, "CA",
                                                                                           " ", "AXS", "A", j + 1, " ", SmoothX[i][j], SmoothY[i][j], SmoothZ[i][j], 1.00, bfactor))

                for j in range(len(SmoothX[i]) - 1):
                    f.write('CONECT %4d %4d\n' % (j + 1, j + 2))

                f.write("TER\n")

            if (write_orig_axis):
                atomstart = 0
                if (write_smooth_axis):
                    atomstart = len(SmoothX[i])

                for j in range(len(RawX[i])):
                    f.write('%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n' % ("ATOM", j + 1 +
                                                                                           atomstart, "O", " ", "AXS", "B", j + 1, " ", RawX[i][j], RawY[i][j], RawZ[i][j], 1.00, 0.00))

                for j in range(len(RawX[i]) - 1):
                    f.write('CONECT %4d %4d\n' %
                            (j + 1 + atomstart, j + 2 + atomstart))

                f.write("TER\n")

            f.write("ENDMDL\n")

        f.close()

    def calculate_curvature_tangent(self, step_range=False, step=None, store_tangent=False):
        """To calculate curvatures and tangent vectors along the helical axis.

        The curvature and tangent vectors are calculated using Frenet-Serret formula.
        The calculated values are stored in ``DNA.data`` dictionary and also in HDF5 file.

        Parameters
        ----------
        step_range : bool
            * ``step_range = True`` : Calculate curvature and tangent vectors for the given range of base-steps
            * ``step_range = False``: Calculate curvature and tangent vectors for entire DNA. If smoothed helical-axis of any base-step will be found to be not available, error will be raised.

        step : list
            List containing lower and higher limit of base-steps range.
                * This option only works with ``step_range=True``.
                * This list should not contain more than two number.
                * First number should be less than second number.

            Example for base-step 4 to 15:
                ``step = [4,15]         # step_range = True``

        store_tangent : bool
            * ``store_tangent = True`` : The calculated tangent vectors will be stored for later use.
            * ``store_tangent = False``:  The calculated tangent vectors will be discarded.

            In case of HDF5 file, calculated tangents will be stroed in this
            file and it will not add cost to memory. However, without HDF5 file,
            storing tangents in ``DNA.data`` will be expansive for memory.

        """

        if not self.smooth_axis:
            raise ValueError(
                "The helical axis is not smooth. At first, smooth the axis using generate_smooth_axis() method as described in http://rjdkmr.github.io/do_x3dna/apidoc.html#dnaMD.DNA.generate_smooth_axis.")

        if (step_range) and (step == None):
            raise ValueError(
                "See, documentation for step  and step_range usage!!!")

        if step_range:
            if (len(step) != 2):
                raise ValueError("See, documentation for step usage!!!")

            if step[0] > step[1]:
                raise ValueError("See, documentation for step usage!!!")

            X, bp_idx = self.get_parameters(
                'helical x-axis smooth', step, bp_range=True)
            Y, bp_idx = self.get_parameters(
                'helical y-axis smooth', step, bp_range=True)
            Z, bp_idx = self.get_parameters(
                'helical z-axis smooth', step, bp_range=True)
        else:
            X, bp_idx = self.get_parameters(
                'helical x-axis smooth', [1, self.num_step], bp_range=True)
            Y, bp_idx = self.get_parameters(
                'helical y-axis smooth', [1, self.num_step], bp_range=True)
            Z, bp_idx = self.get_parameters(
                'helical z-axis smooth', [1, self.num_step], bp_range=True)

        X = np.asarray(X).T
        Y = np.asarray(Y).T
        Z = np.asarray(Z).T

        curvature, tangent = [], []

        for i in range(len(self.time)):

            # Curvature calculation
            xyz = np.vstack((X[i], Y[i], Z[i])).T
            T, N, B, k_temp, t_temp = frenet_serret(xyz)

            curvature.append(k_temp.flatten())

            if(store_tangent):
                tangent.append(T)

        curvature = np.asarray(curvature).T
        for i in range(len(bp_idx)):
            bp_num = str( bp_idx[i]+self.startBP )
            self._set_data(curvature[i], 'bps', bp_num, 'curvature',  scaleoffset=3)

        if(store_tangent):
            tangent = np.asarray(tangent)
            final_tan = []
            for i in range(len(tangent[0])):
                temp = []
                for j in range(len(tangent)):
                    temp.append(tangent[j][i])
                final_tan.append(np.asarray(temp))

            for i in range(len(bp_idx)):
                bp_num = str( bp_idx[i]+self.startBP )
                self._set_data(np.asarray(final_tan[i]), 'bps', bp_num, 'tangent',  scaleoffset=3)

    def calculate_angle_bw_tangents(self, base_step, cumulative=False, masked=False):
        """To calculate angle (Radian) between two tangent vectors of global helical axis.

        Parameters
        ----------
        base_step : 1D list
            List of two base-steps for which angle will be calculated.
            For example: **base_step** = ``[5, 50]`` either of following can be calculated.

                (1) angle between tangent vector 5th and 50th base-steps.
                (2) summation over 44 angles that are formed between adjacent tangent
                    vectors of 5-50 bp DNA segment.

            See below two choose between these two types.

        cumulative : bool)
            ``Default: False``: If it is false, first type of angle is calculated
            otherwise second type of angle is calculated as explained in the above
            example of option ``base_step``.

        masked : bool
            ``Dfault=False``. To skip specific frames/snapshots.
            ``DNA.mask`` array should be set to use this functionality.
            This array contains boolean (either ``True`` or ``False``) value
            for each frame to mask the frames. Presently, mask array is
            automatically generated during :meth:`DNA.generate_smooth_axis` to
            skip those frames where 3D fitting curve was not successfull within
            the given critera.

        Returns
        -------
        angle : 1D array
            Array of calculated angle of length is equal to number of frames.
            When ``masked`` is applied, length of this array can be smaller than
            total number of frames.

        """

        if (len(base_step) != 2):
            raise ValueError("See, documentation for step usage!!!")

        if base_step[0] > base_step[1]:
            raise ValueError("See, documentation for step usage!!!")

        angle = []
        if cumulative:
            angle_tmp = []
            tangent, idx = self.get_parameters(
                'tangent', bp=base_step, bp_range=True, masked=masked)
            for i in range(len(idx) - 1):
                angle_tmp.append(vector_angle(tangent[i], tangent[i + 1]))
            angle = np.asarray(angle_tmp).sum(axis=0)

        else:
            tangent1, idx1 = self.get_parameters(
                'tangent', bp=[base_step[0]], bp_range=False, masked=masked)
            tangent2, idx2 = self.get_parameters(
                'tangent', bp=[base_step[1]], bp_range=False, masked=masked)
            angle = vector_angle(tangent1[0], tangent2[0])

        return np.asarray(angle)

    def calculate_2D_angles_bw_tangents(self, paxis, base_step, masked=True):
        """ To calculate angles (Radian) in 2D plane between two tangent vectors of global helical axis.

        .. list-table::  Principal axis and respective 2D planes
            :widths: 1, 4
            :header-rows: 1
            :name: gui-table

            * - Principal Axis
              - 2D planes

            * - X-axis
              - XY and XZ planes

            * - Y-axis
              - XY and YZ planes

            * - Z-axis
              - XZ and YZ planes

        Parameters
        ----------
        paxis : str
            Principal axis parallel to the DNA axis. For example, IF DNA helical
            axis is parallel or aligned with Z-axis, paxis should be ``'Z'``.
        base_step : 1D list
            List of two base-steps for which angle will be calculated.
            For example: **base_step** = ``[5, 50]`` angle between tangent vector
            5th and 50th base-steps will be calculated.
        masked : bool
            ``Dfault=False``. To skip specific frames/snapshots.
            ``DNA.mask`` array should be set to use this functionality.
            This array contains boolean (either ``True`` or ``False``) value
            for each frame to mask the frames. Presently, mask array is
            automatically generated during :meth:`DNA.generate_smooth_axis` to
            skip those frames where 3D fitting curve was not successfull within
            the given critera.

        Returns
        -------
        angle-one : 1D array
            Array of calculated angle in first plane of length is equal to
            number of frames. When ``masked`` is applied, length of this array
            can be smaller than total number of frames.
        angle-two : 1D array
            Array of calculated angle in second plane of length is equal to
            number of frames. When ``masked`` is applied, length of this array
            can be smaller than total number of frames.

        """

        if (len(base_step) != 2):
            raise ValueError("See, documentation for step usage!!!")

        if base_step[0] > base_step[1]:
            raise ValueError("See, documentation for step usage!!!")

        if paxis == 'Z':
            planes = [0, 1]
        elif paxis == 'X':
            planes = [1, 2]
        elif paxis == 'Y':
            planes = [0, 2]
        else:
            raise ValueError('{0} is not accepted keyword. USE: X Y or Z for DNA-axis.')

        tangent1, idx1 = self.get_parameters('tangent', bp=[base_step[0]], bp_range=False, masked=masked)
        tangent2, idx2 = self.get_parameters('tangent', bp=[base_step[1]], bp_range=False, masked=masked)

        angles = []
        for plane in planes:
            if plane == 0:
                norm = [1, 0, 0]
            elif plane == 1:
                norm = [0, 1, 0]
            else:
                norm = [0, 0, 1]

            temp1 = tangent1[0].T[plane].copy()
            temp2 = tangent2[0].T[plane].copy()

            tangent1[0].T[plane] = np.zeros(tangent1[0].T[plane].shape)
            tangent2[0].T[plane] = np.zeros(tangent2[0].T[plane].shape)

            angle = vector_angle(tangent1[0], tangent2[0], norm=norm)

            tangent1[0].T[plane]  = temp1
            tangent2[0].T[plane]  = temp2

            angles.append(angle)

        return np.asarray(angles[0]), np.asarray(angles[1])

    def calculate_contour_length(self, step_range=False, step=None, masked=True):
        if not self.smooth_axis:
            raise ValueError(
                "The helical axis is not smooth. At first, smooth the axis using generate_smooth_axis() method as described in http://rjdkmr.github.io/do_x3dna/apidoc.html#dnaMD.DNA.generate_smooth_axis.")

        if (step_range) and (step == None):
            raise ValueError(
                "See, documentation for step  and step_range usage!!!")

        if step_range:
            if (len(step) != 2):
                raise ValueError("See, documentation for step usage!!!")

            if step[0] > step[1]:
                raise ValueError("See, documentation for step usage!!!")

            smoothX, bp_idx = self.get_parameters('helical x-axis smooth', step, bp_range=True, masked=masked)
            smoothY, bp_idx = self.get_parameters('helical y-axis smooth', step, bp_range=True, masked=masked)
            smoothZ, bp_idx = self.get_parameters('helical z-axis smooth', step, bp_range=True, masked=masked)
        else:
            smoothX, bp_idx = self.get_parameters('helical x-axis smooth', [2, self.num_step-1], bp_range=True, masked=masked)
            smoothY, bp_idx = self.get_parameters('helical y-axis smooth', [2, self.num_step-1], bp_range=True, masked=masked)
            smoothZ, bp_idx = self.get_parameters('helical z-axis smooth', [2, self.num_step-1], bp_range=True, masked=masked)

        smoothX = np.asarray(smoothX)
        smoothY = np.asarray(smoothY)
        smoothZ = np.asarray(smoothZ)

        distances = []
        for frame in range(len(smoothX[0])):
            pos = np.asarray( [smoothX[:, frame], smoothY[:, frame], smoothZ[:,frame]] )
            distance = np.diagonal( sp_dist.cdist(pos.T, pos.T), offset=1)
            distances.append( distance.sum() )

        return np.asarray(distances)


def checkParametersInputFile(filename):
    """Check the do_x3dna output file and return list of parameters present in the file.
    """
    fin = open(filename, 'r')
    line = fin.readline()
    line2 = fin.readline()
    fin.close()

    temp = re.split('\s+', line)
    temp2 = re.split('\s+', line2)

    if temp[0] == '#Minor':
        return groovesParameters

    if temp[0] == '#Shift':
        return baseStepParameters

    if temp[0] == '#X-disp':
        return helicalBaseStepParameters

    if temp[0] == '#Shear':
        return basePairParameters

    if temp[0] == '#Position':
        return helicalAxisParameters

    if temp2[0] == '#P':
        return helicalRadiusParameters

    if temp2[0] == '#alpha':
        return backboneDihedrals


def setParametersFromFile(dna, filename, parameter, bp=None):
    """Read a specific parameter from the do_x3dna output file.
    It automatically load the input parameter from a file to dna object or HDF5 file.
    It automatically decides from input parameter, what will be format of input file.

    Parameters
    ----------
    dna : :class:`DNA`
        Input :class:`DNA` instance.
    filename : str
        Input filename. This file should be output from do_x3dna.
    parameter : str
        Name of parameter. For details about accepted keywords, see ``parameter`` in the method :meth:`DNA.get_parameters`.
        Note that parameter that are calculated from do_x3dna cannot be used here.
    bp : list
        List containing lower and higher limit of base-pair/step range.
            * This list should not contain more than two number.
            * First number should be less than second number.

        Example for base-pairs/steps 4 to 15:
            ``bp = [4,15]         # step_range = True``

        If ``None``, all base-pairs/steps will be considered.

    """

    gotParamterList = False
    param_type = None
    if isinstance(parameter, list) or isinstance(parameter, np.ndarray):
        gotParamterList = True
        parameter = list(parameter)
        param_type = getParameterType(parameter[0])
    else:
        param_type = getParameterType(parameter)

    if bp is None:
        if param_type == 'bps':
            bp = [dna.startBP, dna.num_step]
        else:
            bp = [dna.startBP, dna.num_bp]

    if len(bp) == 1:
        bp_range = False
    else:
        bp_range = True

    if not gotParamterList:
        tempParamName = parameter
        inputParameter = [ parameter ]
    else:
        tempParamName = parameter[0]
        inputParameter = parameter

    success = False
    if tempParamName in basePairParameters:
        dna.set_base_pair_parameters(filename, bp, parameters=inputParameter, bp_range=bp_range)
        success = True

    if tempParamName in baseStepParameters:
        dna.set_base_step_parameters(filename, bp, parameters=inputParameter, step_range=bp_range, helical=False)
        success = True

    if tempParamName in helicalBaseStepParameters:
        dna.set_base_step_parameters(filename, bp, parameters=inputParameter, step_range=bp_range, helical=True)
        success = True

    if tempParamName in groovesParameters:
        dna.set_major_minor_groove(filename, bp, parameters=inputParameter, step_range=bp_range)
        success = True

    if tempParamName in backboneDihedrals:
        dna.set_backbone_dihedrals(filename, bp, parameters=inputParameter, bp_range=bp_range)
        success = True

    if tempParamName in helicalRadiusParameters:
        dna.set_helical_radius(filename, bp, full=True, bp_range=bp_range)
        success = True

    if tempParamName in helicalAxisParameters:
        if len(bp) == 1:
            raise AssertionError("Axis cannot be read for a single base-step.\n Use a segment spanned over several basepairs.")
        dna.set_helical_axis(filename, step_range=True, step=bp)
        success = True

    if not success:
        raise ValueError ('Not able to set these parameters: {0}... '.format(parameter))


def dev_bps_vs_parameter(dnaRef, bpRef, dnaSubj, bpSubj, parameter, err_type='std', bp_range=True, merge_bp=1, merge_method='mean', tool='g_analyze'):
    """To calculate deviation in the given parameters of a Subject DNA with respect to a Reference DNA along the base-pairs/steps.

    .. note:: Deviation = Reference_DNA(parameter) - Subject_DNA(parameter)

    .. warning:: Number of base-pairs/steps should be similar in reference and subject DNA.

    .. warning::
            To calculate errors by using ``error = 'acf'`` or ``error = 'block'``,
            GROMACS tool ``g_analyze`` or ``gmx analyze`` should be present in ``$PATH``.

    Parameters
    ----------
    dnaRef  : :class:`DNA`
        Reference DNA

    bpRef : 1D list array
        base-pairs or base-steps to consider from Reference DNA
        Example: ::

            bp = [6]                                # bp_range = False
            bp = [4,15]                             # bp_range = True
            bp = range(4,15)                        # bp_range = False
            bp = np.arange(4,15)                    # bp_range = False


    dnaSubj : :class:`DNA`
        Subject DNA. Number of base-pairs in Reference and Subject DNA **should be** same.

    bpSubj : 1D list or array
        base-pairs or base-steps to consider from Reference DNA. Foe more, see above example of ``bpSubj``.

    parameter : str
        Name of a base-pair or base-step or helical base-step parameter
        For details about accepted keywords, see ``parameter`` in the method :meth:`DNA.get_parameters`.

    error_type : str
        Method of error estimation.
        Currently accepted method as follows:

            * ``error = 'std'``   : Standard Deviation
            * ``error = 'acf'``   : Standard error using autocorrelation time (requires: ``g_analyze`` or ``gmx analyze``)
            * ``error = 'block'`` : Standard error using block averaging method (requires: ``g_analyze`` or ``gmx analyze``)

    bp_range : bool
        ``Dfault=True``: As shown above, if ``True``, bp is taken as a range otherwise list or numpy array.

    merge_bp : int
        Number of base-pairs or steps to merge for creating the small DNA segments

    merge_method : str
        Method to calculate the parameter of a DNA segment from local
        parameters of all base-pairs/steps that are between the range
        given through ``bp``.
        Currently accepted keywords are as follows:

            * ``merge_method = mean``: Average of local parameters
            * ``merge_method = sum``: Sum of local parameters

    tool : str
        Gromacs tool ``g_analyze`` or ``gmx analyze`` or ``gmx_mpi analyze`` etc.
        It will be used to calculate autocorrelation time or  block averaging error.
        It should be present in ``$PATH``

    Returns
    -------
    bpRef  : 1D array)
        base-pair/step numbers of reference DNA. If ``merge_bp>1``, middle number will is returned.`
    bpSubj : 1D array
        base-pair/step numbers of subject DNA. If ``merge_bp>1``, middle number will is returned.`
    deviation  : 1D array
        Deviation in the parameter of subject DNA with respect to reference DNA.
    error  : 1D array
        Standard error of respective deviation

    """

    bpRef, RefAvgValue, RefError = dnaRef.get_mean_error(
        bpRef, parameter, err_type=err_type, bp_range=True, merge_bp=merge_bp, merge_method=merge_method, tool=tool)

    bpSubj, SubjAvgValue, SubjError = dnaSubj.get_mean_error(
        bpSubj, parameter, err_type=err_type, bp_range=True, merge_bp=merge_bp, merge_method=merge_method, tool=tool)

    if len(bpRef) != len(bpSubj):
        raise ValueError(
            "Number (%d) of bp/bps/segments in reference DNA does not match with the number (%d) of subject DNA." % (len(bpRef), len(bpSubj)))
        exit(1)

    deviation = RefAvgValue - SubjAvgValue
    error = np.sqrt((RefError**2) + (SubjError**2))

    return bpRef, bpSubj, deviation, error


def dev_parameters_vs_axis(dnaRef, dnaSubj, parameter, bp, axis='Z', bp_range=True, windows=10, err_type='block', tool='g_analyze'):
    """To calculate deviation in the given parameters of a Subject DNA to Reference DNA along the given axis.

    .. note:: Deviation = Reference_DNA(parameter) - Subject_DNA(parameter)

    .. warning::
            To calculate errors by using ``error = 'acf'`` or ``error = 'block'``,
            GROMACS tool ``g_analyze`` or ``gmx analyze`` should be present in ``$PATH``.

    Parameters
    ----------
    dnaRef  : :class:`DNA`
        Reference DNA
    dnaSubj : :class:`DNA`
        Subject DNA. Number of base-pairs in Reference and Subject DNA **should be** same.
    parameter : str
        Name of a base-pair or base-step or helical base-step parameter
        For details about accepted keywords, see ``parameter`` in the method :meth:`DNA.get_parameters`.

    bp : 1D list or array
        base-pairs to analyze
        Example: ::

            bp = [6]                                # bp_range = False
            bp = [4,15]                             # bp_range = True
            bp = range(4,15)                        # bp_range = False
            bp = np.arange(4,15)                    # bp_range = False
            bp = [2,5,6,7,9,12,18]                  # bp_range = False

    bp_range : bool
        ``Dfault=True``: As shown above, if ``True``, bp is taken as a range otherwise list or numpy array.

    axis : str
        Axis along which DNA axis is parallel. Keywords: ``X``, ``Y`` and ``Z``.

    windows  : int
        Number of bins along the axis

    err_type : str
        Method of error estimation.
        Currently accepted method as follows:

            * ``error = 'std'``   : Standard Deviation
            * ``error = 'acf'``   : Standard error using autocorrelation time (requires: ``g_analyze`` or ``gmx analyze``)
            * ``error = 'block'`` : Standard error using block averaging method (requires: ``g_analyze`` or ``gmx analyze``)

    tool : str
        Gromacs tool ``g_analyze`` or ``gmx analyze`` or ``gmx_mpi analyze`` etc.
        It will be used to calculate autocorrelation time or  block averaging error.
        It should be present in ``$PATH``

    Returns
    -------
    deviation : 1D array
        length of no. of windows; Deviation in the parameter for two given DNA
    deviation_error : 1D array
        length of no. of windows; Standard error in deviation fo each window/bin
    axis  : 1D array
        length of no. of windows; average position of window/bin along given axis
    axis_error  : 1D array
        length of no. of windows; Standard error in average position of window/bin along given axis

    """
    RefParam, ref_bp_idx = dnaRef.get_parameters(parameter, bp, bp_range)
    RefAxis, dummy = dnaRef.get_parameters(
        'Helical {0}-axis' .format(axis), bp, bp_range)
    SubjParam, subj_bp_idx = dnaSubj.get_parameters(parameter, bp, bp_range)

    mean_axis = np.mean(RefAxis, axis=1)
    meanRefParam = np.mean(RefParam, axis=1)
    meanSubjParam = np.mean(SubjParam, axis=1)

    maxAxis = np.amax(mean_axis)
    minAxis = np.amin(mean_axis)
    axis_range = (maxAxis - minAxis) / windows

    Ref_param_error = get_error(dnaRef.time, RefParam, len(ref_bp_idx), err_type=err_type, tool=tool)
    Ref_axis_error = get_error(dnaRef.time, RefAxis, len(ref_bp_idx), err_type=err_type, tool=tool)
    subj_param_error = get_error(dnaSubj.time, SubjParam, len(subj_bp_idx), err_type=err_type, tool=tool)

    merged_ref_param = []
    merged_subj_Param = []
    merged_Ref_param_error = []
    merged_Ref_axis_error = []
    merged_subj_param_error = []

    final_axis = []

    for i in range(windows):
        start = minAxis + (i * axis_range)
        end = start + axis_range
        idx = []
        for j in range(len(mean_axis)):
            if((start <= mean_axis[j]) and (end > mean_axis[j])):
                idx.append(j)
        if(len(idx) > 0):
            merged_ref_param.append(meanRefParam[idx])
            merged_subj_Param.append(meanSubjParam[idx])
            final_axis.append(start + (end - start) / 2)
            merged_Ref_param_error.append(Ref_param_error[idx])
            merged_Ref_axis_error.append(Ref_axis_error[idx])
            merged_subj_param_error.append(subj_param_error[idx])

    final_ref_param = []
    final_subj_param = []
    final_ref_param_error = []
    final_ref_axis_error = []
    final_subj_param_error = []
    for i in range(len(merged_ref_param)):
        final_ref_param.append(np.sum(merged_ref_param[i]))
        final_subj_param.append(np.sum(merged_subj_Param[i]))

        final_ref_axis_error.append(
            np.sqrt((merged_Ref_axis_error[i]**2).sum()))
        final_ref_param_error.append(
            np.sqrt((merged_Ref_param_error[i]**2).sum()))
        final_subj_param_error.append(
            np.sqrt((merged_subj_param_error[i]**2).sum()))

    deviation, error = get_deviation(
        final_ref_param, final_ref_param_error, final_subj_param, final_subj_param_error)

    return deviation, error, final_axis, final_ref_axis_error


def get_error(time, x, sets, err_type='block', tool='g_analyze'):
    """To estimate error using block averaging method

    .. warning::
            To calculate errors by using ``error = 'acf'`` or ``error = 'block'``,
            GROMACS tool ``g_analyze`` or ``gmx analyze`` should be present in ``$PATH``.

    Parameters
    ----------
    time : 1D list or array
        :attr:`DNA.time`
    x   : 2D list or array
        Shape of (nset, nframe); where *nset* is number of set and *nframe* is
        total number of frames. *nframe* should be equal to length of time
        list/array
    sets : int
        Number of sets (*nset*)
    err_type : str
        Error estimation by autocorrelation method ``err_type='acf'`` or
        block avearaging method ``err_type='block'``

    Returns
    -------
    error : 1D array
        Of length = number of sets (*nset*)

    """
    for i in range(sets):
        if (len(time) != len(x[i])):
            raise ValueError('\nError: number of frame in time {0} mismatched with {1} of x[{2}]!!\n' .format(
                len(time), len(x[i]), i))

    if not((err_type == 'block') or (err_type == 'acf')):
        print('\nWarning: Method {0} is not implemented. Swtiching to \'acf\'.\n' .format(
            err_type))
        err_type = 'acf'

    error = []
    char_set = string.ascii_lowercase
    name = ''.join(random.sample(string.ascii_lowercase, 10))

    filename = name + '.xvg'
    eefile = 'ee_' + name + '.xvg'
    acfile = 'acf_' + name + '.xvg'

    fout = open(filename, 'w')
    for i in range(len(time)):
        fout.write('{0}' .format(time[i]))
        for j in range(sets):
            fout.write('\t{0}' .format(x[j][i]))
        fout.write("\n")
    fout.close()

    command = '{0} -f {1} -ee {2} -ac {3} -fitfn exp' .format(tool, filename, eefile, acfile)

    p = sub.Popen(command.split(), stdout=sub.PIPE, stderr=sub.PIPE, universal_newlines=True)
    out, outputerror = p.communicate()
    lines = out.split('\n')

    if (err_type == 'block'):
        for line in lines:
            if(re.match('Set', line)):
                temp = line.split()
                error.append(float(temp[3]))

    if (err_type == 'acf'):
        acf_time = []
        for line in lines:
            if(re.match('COR: Correlation time', line)):
                temp = line.split('=')
                acf_time.append(abs(float(temp[1].split()[0])))

        total_time = float(time[-1]) - float(time[0])
        dt = total_time / len(time)
        for i in range(sets):
            if(acf_time[i] >= dt):
                n_indp = total_time / acf_time[i]
                tmp_err = np.std(x[i]) / np.sqrt(n_indp)
            else:
                tmp_err = np.std(x[i]) / np.sqrt(len(time))
            error.append(tmp_err)

    os.remove(filename)
    os.remove(eefile)
    os.remove(acfile)
    if os.path.isfile('fitlog.log'):
        os.remove('fitlog.log')
    return np.array(error)


def get_deviation(Ref, RefErr, x, xerr):
    if (len(Ref) != len(x)):
        print (
            "\nErrori: Number of base pairs/steps mismatch from reference to target!!\n")
        exit(1)
    Ref = np.array(Ref)
    RefErr = np.array(RefErr)
    x = np.array(x)
    xerr = np.array(xerr)

    covariance = np.cov(Ref, x)

    deviation = Ref - x
    error = np.sqrt((RefErr * RefErr) + (xerr * xerr))

    return deviation, error


def frenet_serret(xyz):
    ''' Frenet-Serret Space Curve Invariants

    Taken from "Diffusion Imaging in Python" `DiPy package<http://nipy.org/dipy/>`_

    Calculates the 3 vector and 2 scalar invariants of a space curve defined by vectors r = (x,y,z).  If z is omitted (i.e. the array xyz has shape (N,2), then the curve is only 2D (planar), but the equations are still valid.

    Similar to http://www.mathworks.com/matlabcentral/fileexchange/11169

    In the following equations the prime ($'$) indicates differentiation with respect to the parameter $s$ of a parametrised curve $\mathbf{r}(s)$.

    - $\mathbf{T}=\mathbf{r'}/|\mathbf{r'}|\qquad$ (Tangent vector)}

    - $\mathbf{N}=\mathbf{T'}/|\mathbf{T'}|\qquad$ (Normal vector)

    - $\mathbf{B}=\mathbf{T}\times\mathbf{N}\qquad$ (Binormal vector)

    - $\kappa=|\mathbf{T'}|\qquad$ (Curvature)

    - $\mathrm{\tau}=-\mathbf{B'}\cdot\mathbf{N}$ (Torsion)

    **Arguments:**

            * xyz : array-like shape (N,3)
                            array representing x,y,z of N points in a track

    **Returns:**

                    * T : array shape (N,3)
                            array representing the tangent of the curve xyz
                    * N : array shape (N,3)
                            array representing the normal of the curve xyz
            * B : array shape (N,3)
                    array representing the binormal of the curve xyz
            * k : array shape (N,1)
                    array representing the curvature of the curve xyz
            * t : array shape (N,1)
                    array representing the torsion of the curve xyz

    Examples:

    Create a helix and calculate its tangent, normal, binormal, curvature and torsion

    >>> from dipy.tracking import metrics as tm
    >>> import numpy as np
    >>> theta = 2*np.pi*np.linspace(0,2,100)
    >>> x=np.cos(theta)
    >>> y=np.sin(theta)
    >>> z=theta/(2*np.pi)
    >>> xyz=np.vstack((x,y,z)).T
    >>> T,N,B,k,t=tm.frenet_serret(xyz)

    '''

    def magn(xyz, n=1):
        ''' magnitude of vector


        '''
        mag = np.sum(xyz**2, axis=1)**0.5
        imag = np.where(mag == 0)
        mag[imag] = np.finfo(float).eps

        if n > 1:
            return np.tile(mag, (n, 1)).T

        return mag.reshape(len(mag), 1)

    xyz = np.asarray(xyz)
    n_pts = xyz.shape[0]

    if n_pts == 0:
        raise ValueError('xyz array cannot be empty')

    dxyz = np.gradient(xyz)[0]
    ddxyz = np.gradient(dxyz)[0]

    # Tangent
    T = np.divide(dxyz, magn(dxyz, 3))

    # Derivative of Tangent
    dT = np.gradient(T)[0]

    # Normal
    N = np.divide(dT, magn(dT, 3))

    # Binormal
    B = np.cross(T, N)

    # Curvature
    k = magn(np.cross(dxyz, ddxyz), 1) / (magn(dxyz, 1)**3)

    # Torsion
    #(In matlab was t=dot(-B,N,2))
    t = np.sum(-B * N, axis=1)
    # return T,N,B,k,t,dxyz,ddxyz,dT

    return T, N, B, k, t


def read_data_file(FileName, cols_equal=True):
    infile = open(FileName, 'r')
    data = []
    len_data = 0
    i = 1
    for line in infile:
        line = line.rstrip('\n')
        if not line.strip():
            continue
        if(re.match('#|@', line) == None):
            temp = map(float, line.split())
            if(cols_equal):
                if (i == 1):
                    len_data = len(temp)
                if (len(temp) != len_data):
                    print('WARNING: Number of column mis match at line {0} in {1}; skipping remaining part\n' .format(
                        i, FileName))
                    break
                data.append(temp)
            i = i + 1
    data = np.array(data).T
    return data


def get_idx_of_bp_parameters(bp, parameters, bp_range, startBP=1):
    param_idx = []
    if(bp_range):
        bp_idx = np.arange(bp[0] - startBP, bp[1])
    else:
        bp_idx = np.subtract(bp, startBP)

    if(len(parameters) != 0):
        #param_idx = np.hstack((np.subtract(parameters,1),[parameters[-1]]))
        param_idx = np.subtract(parameters, 1)

    return bp_idx, param_idx


def read_param_file(FileName, parameters, bp, bp_range, word=False, startBP=1):
    """ Read parameters from do_x3dna file.

    It is the main function, which is used to read and extract the parameters
    values from the do_x3dna ouput files.

    Parameters
    ----------
    FileName : str
        Parameter file produced from do_x3dna.
    parameters : list
        List of column indices that has to be extracted. indices here start
        with one.

    bp : 1D list or array
        base-pairs to analyze
        Example: ::

            bp = [6]                                # bp_range = False
            bp = [4,15]                             # bp_range = True
            bp = range(4,15)                        # bp_range = False
            bp = np.arange(4,15)                    # bp_range = False
            bp = [2,5,6,7,9,12,18]                  # bp_range = False
    bp_range : bool
        ``Dfault=True``: As shown above, if ``True``, bp is taken as a range otherwise list or numpy array.
    word : bool
        In some parameters, in place of value, ``'---'`` is present in the file.
        If parameter values contain this, use ``True``.
    startBP : int
        Number ID of first base-pair.

    Returns
    -------
    data : 3D array
        Extracted parameters as a 3D array of shape (bp, parameters, time).
    time : 1D array
        Time of each frame
    """

    sys.stdout.write("\nReading file : %s\n" % FileName)
    sys.stdout.flush()

    def get_frame_data(block, parameters, bp_idx):
        block = np.array(block).T
        temp_data = (block[parameters, :])[:, bp_idx].copy()
        return temp_data

    def get_time(line):
        dummy, temp_time = line.split('=')
        return float( temp_time )

    infile = open(FileName, 'r')
    data = []
    time = []
    frame_number = 0

    bp_idx, param_idx = get_idx_of_bp_parameters(bp, parameters, bp_range, startBP=startBP)

    block = []
    for line in infile:

        # Removing last new line charecter
        line = line.rstrip('\n')

        # Skipping blank/empty line
        if not line.strip():
            continue

        # Getting Time tag and time => Starting of new frame
        if(re.match('# Time', line) != None):

            if((frame_number < 100) and (frame_number % 10 == 0)) or ((frame_number < 1000) and (frame_number % 100 == 0)) or ((frame_number < 10000) and (frame_number % 1000 == 0)) or ((frame_number < 100000) and (frame_number % 10000 == 0)) or ((frame_number < 1000000) and (frame_number % 100000 == 0)):
                sys.stdout.write("\rReading frame %d" % frame_number)
                sys.stdout.flush()

            frame_number += 1

            # if(frame_number==5000):
            # break

            # Getting time
            time.append(get_time(line))

            # Getting parameters/values for base-pairs
            if(len(block) > 0):
                data.append(get_frame_data(block, param_idx, bp_idx))
            block = []
            continue

        # Skipping other lines starting with '#' tag
        if(re.match('#', line) != None):
            continue

        if not word:
            block.append(list(map(float, line.split())))
        else:
            temp = []
            split_line = line.split()
            for word in split_line:
                if word != '---':
                    temp.append(float(word))
                else:
                    temp.append(None)
            block.append(temp)

    # For last frame
    data.append(get_frame_data(block, param_idx, bp_idx))
    block = []
    data_transpose = np.array(data).T

    sys.stdout.write(
        "\nFinishid reading.... Total number of frame read =  %d\n" % frame_number)
    sys.stdout.flush()

    return data_transpose, time


def distance(x, y):
    x = np.asarray(x)
    y = np.asarray(y)
    return np.linalg.norm(x - y)


def vector_angle(x, y, multiple=True, norm=None):
    """Calculate angle between two vector/s

    Parameters
    ----------
    x : 1D or 2D array
        First vector or array of first vectors
    y : 1D or 2D array
        Second vector or array of second vectors
    multiple : bool
        If ``x`` and ``y`` are array of vector then use ``True`` otherwise
        use ``False``
    norm : vector
        Normal vector. useful to determine signed angle.

    Returns
    -------
    value : float or 1D array
        Calculated angle value. If ``x`` and ``y`` are array of vector, then array
        of angle is returned otherwise a single angle value is returned.
    """
    x = np.asarray(x)
    y = np.asarray(y)

    angle = []
    if multiple:
        # Element-wise angle between arrays of vectors. Taken from following links
        # http://codereview.stackexchange.com/questions/54347/calculating-element-wise-the-angles-between-two-lists-of-vectors
        dot_x_y = np.einsum('ij, ij->i', x, y)
        dot_x_x = np.einsum('ij, ij->i', x, x)
        dot_y_y = np.einsum('ij, ij->i', y, y)
        angle = np.arccos(dot_x_y / (np.sqrt(dot_x_x) * np.sqrt(dot_y_y)))

        if norm is not None:
            cross_x_y = np.cross(x, y)
            sign = np.einsum('ij, ij->i', cross_x_y, np.tile(norm, (x.shape[0], 1)))
            angle = angle * np.sign(sign)
    else:
        # Angle between two vectors
        dot = np.dot(x, y)
        cross = np.cross(x, y)
        cross_modulus = np.sqrt((cross * cross).sum())
        angle = np.arctan2(cross_modulus, dot)
    return angle

####################################


def fit_axis(bp_idx, nframe, RawX, RawY, RawZ, smooth, spline, fill_point, cut_off_angle):

    orig_x = RawX.copy()
    orig_y = RawY.copy()
    orig_z = RawZ.copy()

    orig_s = smooth

    xsmooth, ysmooth, zsmooth = [], [], []

    bsmooth = False
    del_idx = []
    count = 0
    scount = 0
    mask = False
    while not bsmooth:

        count += 1

        if count > 4:
            mask = False
            sys.stdout.write('\n|frame:{0:>10}| WARNING: Fitting failed with \"smooth = {1}\"; Trying with \"smooth = {2}\".....\n' .format(
                nframe, smooth, smooth + 100))
            sys.stdout.flush()

            smooth = smooth + 100
            count = 1

            del orig_x
            del orig_y
            del orig_z

            orig_x = RawX.copy()
            orig_y = RawY.copy()
            orig_z = RawZ.copy()

        points = fill_point * len(orig_x)

        nest = -1
        tckp, u = splprep(
            [orig_x, orig_y, orig_z], s=smooth, k=spline, nest=nest)

        xnew, ynew, znew = splev(np.linspace(0, 1, points), tckp)

        new_axis = np.array([xnew, ynew, znew]).T

        angle = []
        dist = []
        del_idx = []
        last_idx = len(orig_x) - 1

        for nbp in range(len(bp_idx)):
            start = nbp * fill_point
            end = start + fill_point
            xsmooth.append(xnew[start:end].mean())
            ysmooth.append(ynew[start:end].mean())
            zsmooth.append(znew[start:end].mean())

        tmp_vec1, tmp_vec2 = [], []
        for j in range(1, len(xsmooth) - 1):
            prev = np.array([xsmooth[j - 1], ysmooth[j - 1], zsmooth[j - 1]])
            curr = np.array([xsmooth[j], ysmooth[j], zsmooth[j]])
            nex = np.array([xsmooth[j + 1], ysmooth[j + 1], zsmooth[j + 1]])
            tmp_vec1.append(prev - curr)
            tmp_vec2.append(curr - nex)
        angle = np.degrees(vector_angle(tmp_vec1, tmp_vec2))
        del tmp_vec1
        del tmp_vec2

        for j in range(1, len(orig_x) - 1):
            prev = np.array([orig_x[j - 1], orig_y[j - 1], orig_z[j - 1]])
            curr = np.array([orig_x[j], orig_y[j], orig_z[j]])
            nex = np.array([orig_x[j + 1], orig_y[j + 1], orig_z[j + 1]])
            dist.append(distance(prev, curr) + distance(curr, nex))

        for j in range(len(angle)):

            if angle[j] > cut_off_angle and not angle[j] > (180 - cut_off_angle):

                del orig_x
                del orig_y
                del orig_z
                bsmooth = False

                max_idx = np.argsort(dist)[::-1]

                for k in range(count):
                    del_idx.append(max_idx[k] + 1)
                    del_idx.append(max_idx[k] + 2)
                    if max_idx[k] == 0:
                        del_idx.append(0)
                    if max_idx[k] == last_idx:
                        del_idx.append(last_idx)

                del_idx = list(set(del_idx))

                sys.stdout.write('\r|frame:{0:>10}| WARNING: Bending angle [{1}-{2}-{3}] = {4:.2f} is more than cut-off angle {5};\n                     Four maximum distances between three adjacent axis positions = ({6[0]:.1f}, {6[1]:.1f}, {6[2]:.1f}, {6[3]:.1f});\n                     Deleting {7} original helical axis positions to remove possible fitting artifact...\n' .format(
                    nframe, j - 1, j, j + 1, angle[j], cut_off_angle, np.sort(dist)[::-1], del_idx))
                sys.stdout.flush()
                mask = True

                if smooth >= 10000:
                    sys.stdout.write('\n\n|frame:{0:>10}| WARNING: Maximum Bending Angle = {1:.2f} at index {2}, which might be artifect. Please, check by visualizing PDB trajectory file...\n\n' .format(
                        nframe, angle[j], j))
                    sys.stdout.flush()
                    bsmooth = True
                    mask = True
                    break

                orig_x = np.delete(RawX.copy(), del_idx)
                orig_y = np.delete(RawY.copy(), del_idx)
                orig_z = np.delete(RawZ.copy(), del_idx)

                xsmooth, ysmooth, zsmooth = [], [], []

                break
            else:
                bsmooth = True

        del angle
        del dist

    if orig_s != smooth:
        sys.stdout.write('\n')

    return xsmooth, ysmooth, zsmooth, mask
#################################################
