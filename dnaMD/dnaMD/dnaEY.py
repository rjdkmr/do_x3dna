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

import numpy as np

from . import dnaMD


class dnaEY:
    """ Elastic Properties Class to calculate elasticity and deformation free energy

    **To initialize this class:** ::

        dna = dnaEY(60, 'BST', filename='dna.h5')       # 60 is the number of basepairs

    This class also contains several methods (functions) that are discussed in following sections.

    Parameters
    ---------
    num_bp : int
        Number of base-pairs in DNA
    esType : str
        Elastic Properties type: ``'BST'`` or ``'ST'``

        (1) ``'BST'``: Bending-Stretching-Twisting --- All motions are considered
        (2) ``'ST'`` : Stretching-Twisting --- Bending motions are ignored.

        .. warning:: For accurate calculation of bending motions, DNA structures in trajectory must be superimposed
                on a reference structure (**See Publication's Method Section**).

    filename : str
        HDF5 (.h5) file containing DNA data.
        In case of ``esType='BST'``, it must contains coordinate of global helical axis (done by
        :meth:`dnaMD.generate_smooth_axis` or 'dnaMD axisCurv').

    startBP : int
        Number ID of first basepair.


    Attributes
    ----------
    num_bp : int
        Number of base-pairs in DNA
    esType : str
        Elastic Properties type: ``'BST'`` or ``'ST'``

        (1) ``'BST'``: Bending-Stretching-Twisting --- All motions are considered
        (2) ``'ST'`` : Stretching-Twisting --- Bending motions are ignored.

        .. warning:: For accurate calculation of bending motions, DNA structures in trajectory must be superimposed
                on a reference structure (**See Publication's Method Section**).

    filename : str
        HDF5 (.h5) file containing DNA data.
        In case of ``esType='BST'``, it must contains coordinate of global helical axis (done by
        :meth:`dnaMD.generate_smooth_axis` or ``dnaMD axisCurv``).

    startBP : int
        Number ID of first basepair.

    esMatrix : dict
        Dictionary of Elastic Matrices. When a combination of new DNA segment and Trajectory segment is given,
        its Elastic Matrix is calculated and stored in this dictionary. It reduces time, when same DNA segment and
        Trajectory segment is given to calculate elastic properties.

        .. note:: This and :attr:`dnaEY.minimumPoint` dictionary can be saved as a json file for later use.

    minimumPoint : dict
        Dictionary of minimum points on energy landscape. The dictionary keywords are similar to that of
        :attr:`dnaEY.esMatrix`. It is generated to speed up the calculation.
        It contains minimum (average) points for the motions. It is used in the calculation of deformation free energy.

        .. note:: This and :attr:`dnaEY.esMatrix` dictionary can be saved as a json file for later use.

    """

    def __init__(self, num_bp, esType='ST', filename=None, startBP=1):
        """Initialize Elastic Properties class

        """
        self.dna = dnaMD.DNA(num_bp, filename=filename, startBP=startBP)

        if esType not in ['BST', 'ST']:
            raise ValueError('Accepted keywords are BST and ST !!!')
        self.esType = esType

        self.esMatrix = dict()             # bp[0]-bp[1]-frame[0]-frame[1]
        self.minimumPoint = dict()         # bp[0]-bp[1]-frame[0]-frame[1]

        self.enGlobalTypes = ['all', 'diag', 'b1', 'b2', 'stretch', 'twist', 'st_coupling',
                        'bend', 'bs_coupling', 'bt_coupling', 'bb_coupling', 'st', 'bb', 'bs', 'bt', ]

        self.enLocalTypes = ['all', 'diag', 'shift', 'slide', 'rise', 'tilt', 'roll', 'twist',
                             'x-disp', 'y-disp', 'h-rise', 'inclination', 'tip', 'h-twist' ]


    def getElasticMatrix(self, inputArray):
        inputArray = np.array( inputArray )

        # Calculation of covariance matrix
        CovMat = np.cov(inputArray, bias=True)

        # Change to a matrix object
        CovMat = np.matrix(CovMat)

        # Inverse of the covariance matrix
        InvCovMat = CovMat.I

        return np.asarray(InvCovMat)

    def _validateFrames(self, frames):

        if frames is None:
            frames = [0, -1]
        else:
            if (len(frames) != 2):
                raise ValueError("frames should be a list containing lower and higher limit. See, documentation!!!")

            if frames[1] != -1 and frames[0] > frames[1]:
                raise ValueError("frames should be a list containing lower and higher limit. See, documentation!!!")

        return frames

    def _extractGlobalParameters(self, dna, bp, frames=None, paxis='Z', masked=False):
        """Extract the parameters for calculations


        .. currentmodule:: dnaMD

        Parameters
        ----------
        dna : :class:`dnaMD.DNA`
            Input :class:`dnaMD.DNA` instance.

        bp : list
            List of two base-steps forming the DNA segment.
            For example: with ``bp=[5, 50]``, 5-50 base-step segment will be considered.

        frames : list
            List of two trajectory frames between which parameters will be extracted. It can be used to select portions
            of the trajectory. For example, with ``frames=[100, 1000]``, 100th to 1000th frame of the trajectory will be
            considered.

        paxis : str
            Axis parallel to global helical-axis(``'X'``, or ``'Y'`` or ``'Z'``). Only require when bending motions are
            included in the calculation.

        masked : bool
            ``Default=False``. To skip specific frames/snapshots.
            ``DNA.mask`` array should be set to use this functionality.
            This array contains boolean (either ``True`` or ``False``) value
            for each frame to mask the frames. Presently, mask array is
            automatically generated during :meth:`dnaMD.DNA.generate_smooth_axis` to
            skip those frames where 3D fitting curve was not successful within
            the given criteria.


        Returns
        -------
        time : numpy.ndarray
            1D numpy array of shape (nframes) containing time
        array : numpy.ndarray
            2D numpy array of shape (parameters count, nframes) containing extracted parameters.
        """

        frames = self._validateFrames(frames)

        if (len(bp) != 2):
            raise ValueError("bp should be a list containing first and last bp of a segment. See, documentation!!!")

        if bp[0] > bp[1]:
            raise ValueError("bp should be a list containing first and last bp of a segment. See, documentation!!!")

        time, clen = dna.time_vs_parameter('h-rise', bp=bp, merge=True, merge_method='sum', masked=masked)
        clen = np.asarray(clen) * 0.1  # conversion to nm

        time, htwist = dna.time_vs_parameter('h-twist', bp=bp, merge=True, merge_method='sum', masked=masked)
        htwist = np.deg2rad(htwist)  # Conversion to radian

        angleOne, angleTwo = None, None
        if self.esType=='BST':
            angleOne, angleTwo = dna.calculate_2D_angles_bw_tangents(paxis, bp, masked=masked)

            # Rarely there are nan during angle calculation, remove those nan
            nanInOne = np.isnan(angleOne[frames[0]:frames[1]])
            nanInTwo = np.isnan(angleTwo[frames[0]:frames[1]])
            notNan = ~(nanInOne + nanInTwo)
            notNanIdx = np.nonzero(notNan)

            if frames[1] == -1:
                array = np.array([angleOne[frames[0]:][notNanIdx], angleTwo[frames[0]:][notNanIdx],
                                  clen[frames[0]:][notNanIdx], htwist[frames[0]:][notNanIdx]])
                time = (time[frames[0]:])[notNanIdx]
            else:

                array = np.array([angleOne[frames[0]:frames[1]][notNanIdx], angleTwo[frames[0]:frames[1]][notNanIdx],
                                  clen[frames[0]:frames[1]][notNanIdx], htwist[frames[0]:frames[1]][notNanIdx]])
                time = (time[frames[0]:frames[1]])[notNanIdx]

        else:
            if frames[1] == -1:
                array = np.array([clen[frames[0]:], htwist[frames[0]:]])
                time = time[frames[0]:]
            else:
                array = np.array([clen[frames[0]:frames[1]], htwist[frames[0]:frames[1]]])
                time = time[frames[0]:frames[1]]

        return time, array

    def _extractLocalParameters(self, dna, bp, helical=False, frames=None):
        """Extract the local parameters for calculations


        .. currentmodule:: dnaMD

        Parameters
        ----------
        dna : :class:`dnaMD.DNA`
            Input :class:`dnaMD.DNA` instance.

        bp : list
            List of two base-steps forming the DNA segment.
            For example: with ``bp=[5, 50]``, 5-50 base-step segment will be considered.

        frames : list
            List of two trajectory frames between which parameters will be extracted. It can be used to select portions
            of the trajectory. For example, with ``frames=[100, 1000]``, 100th to 1000th frame of the trajectory will be
            considered.

        helical : bool
            If ``helical=True``, helical base-step parameters are extracted. Otherwise, by default, base-step parameters
            are extracted.

        Returns
        -------
        time : numpy.ndarray
            1D numpy array of shape (nframes) containing time
        array : numpy.ndarray
            2D numpy array of shape (6, nframes) containing extracted parameters.
        """

        frames = self._validateFrames(frames)

        if (len(bp) != 2):
            raise ValueError("bp should be a list containing first and last bp of a segment. See, documentation!!!")

        if bp[0] > bp[1]:
            raise ValueError("bp should be a list containing first and last bp of a segment. See, documentation!!!")

        if (bp[1] - bp[0]) > 4:
            print("WARNING: this is a local property and therefore, longer than 4 base-step may not be suitable..." )

        if not helical:
            time, shift = dna.time_vs_parameter('shift', bp=bp, merge=True, merge_method='sum')
            shift = np.asarray(shift) * 0.1  # conversion to nm

            time, slide = dna.time_vs_parameter('slide', bp=bp, merge=True, merge_method='sum')
            slide = np.asarray(slide) * 0.1  # conversion to nm

            time, rise = dna.time_vs_parameter('rise', bp=bp, merge=True, merge_method='sum')
            rise = np.asarray(rise) * 0.1  # conversion to nm

            time, tilt = dna.time_vs_parameter('tilt', bp=bp, merge=True, merge_method='sum')

            time, roll = dna.time_vs_parameter('roll', bp=bp, merge=True, merge_method='sum')

            time, twist = dna.time_vs_parameter('twist', bp=bp, merge=True, merge_method='sum')

            if frames[1] == -1:
                array = np.array([shift[frames[0]:], slide[frames[0]:], rise[frames[0]:],
                                  tilt[frames[0]:], roll[frames[0]:], twist[frames[0]:]])
                time = time[frames[0]:]
            else:
                array = np.array([shift[frames[0]:frames[1]], slide[frames[0]:frames[1]], rise[frames[0]:frames[1]],
                                  tilt[frames[0]:frames[1]], roll[frames[0]:frames[1]], twist[frames[0]:frames[1]]])
                time = time[frames[0]:frames[1]]

        else:
            time, x_disp = dna.time_vs_parameter('x-disp', bp=bp, merge=True, merge_method='sum')
            x_disp = np.asarray(x_disp) * 0.1  # conversion to nm

            time, y_disp = dna.time_vs_parameter('y-disp', bp=bp, merge=True, merge_method='sum')
            y_disp = np.asarray(y_disp) * 0.1  # conversion to nm

            time, h_rise = dna.time_vs_parameter('h-rise', bp=bp, merge=True, merge_method='sum')
            h_rise = np.asarray(h_rise) * 0.1  # conversion to nm

            time, inclination = dna.time_vs_parameter('inclination', bp=bp, merge=True, merge_method='sum')

            time, tip = dna.time_vs_parameter('tip', bp=bp, merge=True, merge_method='sum')

            time, h_twist = dna.time_vs_parameter('h-twist', bp=bp, merge=True, merge_method='sum')

            if frames[1] == -1:
                array = np.array([x_disp[frames[0]:], y_disp[frames[0]:], h_rise[frames[0]:],
                                  inclination[frames[0]:], tip[frames[0]:], h_twist[frames[0]:]])
                time = time[frames[0]:]
            else:
                array = np.array([x_disp[frames[0]:frames[1]], y_disp[frames[0]:frames[1]], h_rise[frames[0]:frames[1]],
                                  inclination[frames[0]:frames[1]], tip[frames[0]:frames[1]], h_twist[frames[0]:frames[1]]])
                time = time[frames[0]:frames[1]]

        return time, array

    def getStretchTwistBendModulus(self, bp, frames=None, paxis='Z', masked=True, matrix=False):
        r"""Calculate Bending-Stretching-Twisting matrix

        It calculate elastic matrix and modulus matrix.

        .. math::

            \text{modulus matrix} = 4.1419464 \times \begin{bmatrix}
            K_{Bx}       & K_{Bx,By} & K_{Bx,S} & K_{Bx,T} \\
            K_{Bx,By}    & K_{By}    & K_{By,S} & K_{By,T} \\
            K_{Bx,S}     & K_{By,S}  & K_{S}    & K_{S,T} \\
            K_{Bx,T}     & K_{Bx,T}  & K_{S,T}  & K_{T}
            \end{bmatrix} \times L_0


        .. currentmodule:: dnaMD

        Parameters
        ----------
        bp : list
            List of two base-steps forming the DNA segment.
            For example: with ``bp=[5, 50]``, 5-50 base-step segment will be considered.

        frames : list
            List of two trajectory frames between which parameters will be extracted. It can be used to select portions
            of the trajectory. For example, with ``frames=[100, 1000]``, 100th to 1000th frame of the trajectory will be
            considered.

        paxis : str
            Axis parallel to global helical-axis(``'X'``, or ``'Y'`` or ``'Z'``). Only require when bending motions are
            included in the calculation.

        masked : bool
            ``Default=False``. To skip specific frames/snapshots.
            ``DNA.mask`` array should be set to use this functionality.
            This array contains boolean (either ``True`` or ``False``) value
            for each frame to mask the frames. Presently, mask array is
            automatically generated during :meth:`dnaMD.DNA.generate_smooth_axis` to
            skip those frames where 3D fitting curve was not successful within
            the given criteria.

        matrix : bool
            If it is ``True``, elastic constant matrix will be returned. Otherwise, by default modulus matrix will be
            returned.

        Return
        ------
        mean : numpy.ndarray
            Value of bending angles, contour length and twist angle (as 1D array) at which energy is zero. Minimum point
            on free energy landscape.

            .. math::
                \begin{bmatrix}
                    \theta^{x}_0       & \theta^{y}_0 & L_0 & \phi_0
                \end{bmatrix}

        result : numpy.ndarray
            Either elastic matrix or modulus matrix depending on ``matrix`` value.
        """


        if self.esType == 'ST':
            raise KeyError(' Use dnaEY.getStretchTwistModulus for Stretching-Twisting modulus.')

        frames = self._validateFrames(frames)

        name = '{0}-{1}-{2}-{3}'.format(bp[0], bp[1], frames[0], frames[1])

        if name not in self.esMatrix:
            time, array = self._extractGlobalParameters(self.dna, bp, frames=frames, paxis=paxis, masked=masked)
            mean = np.mean(array, axis=1)
            esMatrix = np.asarray(self.getElasticMatrix(array))
            self.esMatrix[name] = esMatrix
            self.minimumPoint[name] = mean
        else:
            esMatrix = self.esMatrix[name]
            mean = self.minimumPoint[name]

        if not matrix:
            result = 4.1419464 * np.array(esMatrix) * mean[2]       # Calculate modulus
        else:
            result = esMatrix

        return mean, result

    def getStretchTwistModulus(self, bp, frames=None, masked=False, matrix=False):
        r"""Calculate Stretching-Twisting matrix

        It calculate elastic matrix and modulus matrix.

        .. math::
            \text{modulus matrix} = 4.1419464 \times \begin{bmatrix}
            K_{S}    & K_{S,T} \\
            K_{S,T}  & K_{T}
            \end{bmatrix} \times L_0


        .. currentmodule:: dnaMD

        Parameters
        ----------
        bp : list
            List of two base-steps forming the DNA segment.
            For example: with ``bp=[5, 50]``, 5-50 base-step segment will be considered.

        frames : list
            List of two trajectory frames between which parameters will be extracted. It can be used to select portions
            of the trajectory. For example, with ``frames=[100, 1000]``, 100th to 1000th frame of the trajectory will be
            considered.

        masked : bool
            ``Default=False``. To skip specific frames/snapshots.
            ``DNA.mask`` array should be set to use this functionality.
            This array contains boolean (either ``True`` or ``False``) value
            for each frame to mask the frames. Presently, mask array is
            automatically generated during :meth:`dnaMD.DNA.generate_smooth_axis` to
            skip those frames where 3D fitting curve was not successful within
            the given criteria.

        matrix : bool
            If it is ``True``, elastic constant matrix will be returned. Otherwise, by default modulus matrix will be
            returned.

        Return
        ------
        mean : numpy.ndarray
            Value of bending angles, contour length and twist angle (as 1D array) at which energy is zero. Minimum point
            on free energy landscape.

            .. math::
                \begin{bmatrix}
                    L_0 & \phi_0
                \end{bmatrix}

        result : numpy.ndarray
            Either elastic matrix or modulus matrix depending on ``matrix`` value.
        """

        if self.esType == 'BST':
            raise KeyError(' Use dnaEY.getStretchTwistBendModulus for Bending-Stretching-Twisting modulus.')

        frames = self._validateFrames(frames)

        name = '{0}-{1}-{2}-{3}'.format(bp[0], bp[1], frames[0], frames[1])

        if name not in self.esMatrix:
            time, array = self._extractGlobalParameters(self.dna, bp, frames=frames, masked=masked)
            mean = np.mean(array, axis = 1)
            esMatrix = self.getElasticMatrix(array)
            self.esMatrix[name] = esMatrix
            self.minimumPoint[name] = mean
        else:
            esMatrix = self.esMatrix[name]
            mean = self.minimumPoint[name]

        if not matrix:
            result = 4.1419464 * np.array(esMatrix) * mean[0]       # Calculate modulus
        else:
            result = esMatrix

        return mean, result

    def getModulusByTime(self, bp, frameGap, masked=False, paxis='Z'):
        r"""Calculate moduli as a function of time for convergence check

        It can be used to obtained elastic moduli as a function of time to check their convergence.

        .. note:: Elastic properties cannot be calculated using a single frame because fluctuations are required.
                  Therefore, here time means trajectory between zero time to given time.

        When ``esType='BST'``, following is obtained:

            1) bending-1

            2) bending-2

            3) Stretching

            4) Twisting

            5) bending-1-bending-2

            6) bending-2-stretching

            7) Stretching-Twisting

            8) bending-1-stretching

            9) bending-2-Twisting

            10) bending-1-twisting


        When ``esType='ST'``, following is obtained:

            1) Stretching

            2) Twisting

            3) Stretching-Twisting

        .. currentmodule:: dnaMD

        Parameters
        ----------
        bp : list
            List of two base-steps forming the DNA segment.
            For example: with ``bp=[5, 50]``, 5-50 base-step segment will be considered.

        frameGap : int
            How many frames to skip for next calculation. this option will determine the
            time-gap between each calculation. Lower the number, slower will be the calculation.

        masked : bool
            ``Default=False``. To skip specific frames/snapshots.
            ``DNA.mask`` array should be set to use this functionality.
            This array contains boolean (either ``True`` or ``False``) value
            for each frame to mask the frames. Presently, mask array is
            automatically generated during :meth:`dnaMD.DNA.generate_smooth_axis` to
            skip those frames where 3D fitting curve was not successful within
            the given criteria.

        paxis : str
            Axis parallel to global helical-axis(``'X'``, or ``'Y'`` or ``'Z'``). Only require when bending motions are
            included in the calculation.


        Returns
        -------
        time : numpy.ndarray
            1D array containing time values.

        modulus : numpy.ndarray
            2D array of shape (10, nframes) or (3, nframes) depending on the ``esType``. Order of modulus in this array
            is same as to the above listed order.

        """


        length = len(self.dna.time[:])

        time, modulus = [], []
        for i in range(frameGap, length, frameGap):

            props = None
            if self.esType == 'BST':
                mean, modulus_t = self.getStretchTwistBendModulus(bp, frames=[0, i], paxis=paxis, masked=True)

                props = np.diagonal(modulus_t, offset=0)
                props = np.hstack((props, np.diagonal(modulus_t, offset=1)))
                props = np.hstack((props, np.diagonal(modulus_t, offset=2)))
                props = np.hstack((props, np.diagonal(modulus_t, offset=3)))

            if self.esType == 'ST':
                mean, modulus_t = self.getStretchTwistModulus(bp, frames=[0, i], masked=masked)

                props = np.diagonal(modulus_t, offset=0)
                props = np.hstack((props, np.diagonal(modulus_t, offset=1)))

            time.append(self.dna.time[i])
            modulus.append(props)

        modulus = np.asarray(modulus)

        return time, modulus.T

    def getGlobalDeformationEnergy(self, bp, complexDna, freeDnaFrames=None, boundDnaFrames=None, paxis='Z', which='all', masked=False):
        r"""Deformation energy of the input DNA using Global elastic properties

        It can be used to calculated deformation energy of a input DNA with reference to the DNA present in the current
        object.

        The deformation free energy is calculated using elastic matrix as follows

        .. math::

            G = \frac{1}{2L_0}\mathbf{xKx^T}


        .. math::
            \mathbf{x} =  \begin{bmatrix}
                    (\theta^{x} - \theta^{x}_0)    & (\theta^{y} - \theta^{y}_0) & (L - L_0) & (\phi - \phi_0)
                          \end{bmatrix}


        .. currentmodule:: dnaMD


        Parameters
        ----------
        bp : list
            List of two base-steps forming the DNA segment.
            For example: with ``bp=[5, 50]``, 5-50 base-step segment will be considered.
        complexDna : :class:`dnaMD.DNA`
            Input :class:`dnaMD.DNA` instance for which deformation energy will be calculated.
        freeDnaFrames : list
            To select a trajectory segment of current (free) DNA data.
            List of two trajectory frames between which parameters will be extracted. It can be used to select portions
            of the trajectory. For example, with ``frames=[100, 1000]``, 100th to 1000th frame of the trajectory will be
            considered.
        boundDnaFrames : list
            To select a trajectory segment of input (bound) DNA data.
            List of two trajectory frames between which parameters will be extracted. It can be used to select portions
            of the trajectory. For example, with ``frames=[100, 1000]``, 100th to 1000th frame of the trajectory will be
            considered.
        paxis : str
            Axis parallel to global helical-axis(``'X'``, or ``'Y'`` or ``'Z'``). Only require when bending motions are
            included in the calculation.
        which : all
            For which motions, energy should be calculated.

            Following keywords are available:
                * ``'all'`` : Use entire elastic matrix -- all motions with their coupling
                * ``'diag'`` : Use diagonal of elastic matrix -- all motions but no coupling
                * ``'b1'`` : Only bending-1 motion
                * ``'b2'`` : Only bending-2 motion
                * ``'stretch'`` : Only stretching motion
                * ``'twist'`` : Only Twisting motions
                * ``'st_coupling'`` : Only stretch-twist coupling motion
                * ``'bs_coupling'`` : Only Bending and stretching coupling
                * ``'bt_coupling'`` : Only Bending and Twisting coupling
                * ``'bb_coupling'`` : Only bending-1 and bending-2 coupling
                * ``'bend'`` : Both bending motions with their coupling
                * ``'st'`` : Stretching and twisting motions with their coupling
                * ``'bs'`` : Bending (b1, b2) and stretching motions with their coupling
                * ``'bt'`` : Bending (b1, b2) and twisting motions with their coupling
        masked : bool
            ``Default=False``. To skip specific frames/snapshots.
            ``DNA.mask`` array should be set to use this functionality.
            This array contains boolean (either ``True`` or ``False``) value
            for each frame to mask the frames. Presently, mask array is
            automatically generated during :meth:`dnaMD.DNA.generate_smooth_axis` to
            skip those frames where 3D fitting curve was not successful within
            the given criteria.


        Returns
        -------
        time : numpy.ndarray
            1D array containing time values.

        energy : numpy.ndarray
            1D array of shape (nframes) containing energy value for query DNA.

        """
        if which not in self.enGlobalTypes:
            raise ValueError('{0} is not a supported energy keywords.\n Use any of the following: \n {1}'.format(
                which, self.enGlobalTypes))

        if self.esType == 'BST':
            means, esMatrix = self.getStretchTwistBendModulus(bp, frames=freeDnaFrames, masked=masked, matrix=True, paxis=paxis)
        else:
            means, esMatrix = self.getStretchTwistModulus(bp, frames=freeDnaFrames, masked=masked, matrix=True)

        esMatrix = 2.5 * esMatrix          # Convert kT to kJ/mol
        time, array = self._extractGlobalParameters(complexDna, bp, frames=boundDnaFrames, paxis=paxis, masked=masked)

        energy = []
        diffs = []
        for i in range(array[0].shape[0]):
            vec = array[:, i]
            diff = vec - means
            diffs.append(diff)
            if self.esType == 'BST':
                t_energy = self._calcEnergyBendStretchTwist(diff, esMatrix, which)
            else:
                t_energy = self._calcEnergyStretchTwist(diff, esMatrix, which)

            energy.append(t_energy)

        return time, np.asarray(energy)

    def _calcEnergyStretchTwist(self, diff, es, which):
        r"""Calculate energy for ``estype='ST'`` using a difference vector.

        It is called in :meth:`dnaEY.getGlobalDeformationEnergy` for energy calculation of each frame.

        Parameters
        ----------
        diff : numpy.ndarray
            Array of difference between minimum and current parameter values.

            .. math::

                \mathbf{x} =  \begin{bmatrix}
                    (L_i - L_0) & (\phi_i - \phi_0)
                          \end{bmatrix}

        es : numpy.ndarray
            Elastic matrix. See in :meth:`dnaEY.getStretchTwistModulus` about elastic matrix.

        which : str
            For which type of motions, energy will be calculated.
            See ``which`` parameter in :meth:`dnaEY.getGlobalDeformationEnergy` for keywords.

        Return
        ------
        energy : float
            Deformation free energy value

        """

        if which not in self.enGlobalTypes[:5]:
            raise ValueError('{0} is not a supported energy keywords.\n Use any of the following: \n {1}'.format(
                which, self.enGlobalTypes))

        energy = None

        if which == 'all':
            temp = np.matrix(diff)
            energy = 0.5 * ((temp * es) * temp.T)
            energy = energy[0,0]

        if which == 'diag':
            energy = 0.5 * ((diff[0] ** 2 * es[0][0])
                            + (diff[1] ** 2 * es[1][1]))

        if which == 'stretch':
            energy = 0.5 * (diff[0] ** 2 * es[0][0])

        if which == 'twist':
            energy = 0.5 * (diff[1] ** 2 * es[1][1])

        if which == 'st_coupling':
            energy = 0.5 * (diff[0] * diff[1] * es[0][1])

        return energy

    def _calcEnergyBendStretchTwist(self, diff, es, which):
        r"""Calculate energy for ``esType='BST'`` using a difference vector.

        It is called in :meth:`dnaEY.getGlobalDeformationEnergy` for energy calculation of each frame.

        Parameters
        ----------
        diff : numpy.ndarray
            Array of difference between minimum and current parameter values.

            .. math::

                \mathbf{x} =  \begin{bmatrix}
                    (\theta^{x}_{i} - \theta^{x}_0)    & (\theta^{y}_{i} - \theta^{y}_0) & (L_i - L_0) & (\phi_i - \phi_0)
                          \end{bmatrix}

        es : numpy.ndarray
            Elastic matrix. See in :meth:`dnaEY.getStretchTwistBendModulus` about elastic matrix.

        which : str
            For which type of motions, energy will be calculated.
            see ``which`` parameter in :meth:`dnaEY.getGlobalDeformationEnergy` for keywords.

        Return
        ------
        energy : float
            Deformation free energy value

        """

        if which not in self.enGlobalTypes:
            raise ValueError('{0} is not a supported energy keywords.\n Use any of the following: \n {1}'.format(
                which, self.enGlobalTypes))

        energy = None

        if which == 'all':
            temp = np.matrix(diff)
            energy = 0.5 * ((temp * es) * temp.T)
            energy = energy[0,0]

        if which == 'diag':
            energy = 0.5 * ((diff[0] ** 2 * es[0][0])
                            + (diff[1] ** 2 * es[1][1])
                            + (diff[2] ** 2 * es[2][2])
                            + (diff[3] ** 2 * es[3][3]))

        if which == 'bend':
            energy = 0.5 * ((diff[0] ** 2 * es[0][0])
                            + (diff[1] ** 2 * es[1][1])
                            + (diff[0] * diff[1] * es[0][1]))

        if which == 'b1':
            energy = 0.5 * (diff[0] ** 2 * es[0][0])

        if which == 'b2':
            energy = 0.5 * (diff[1] ** 2 * es[1][1])

        if which == 'stretch':
            energy = 0.5 * (diff[2] ** 2 * es[2][2])

        if which == 'twist':
            energy = 0.5 * (diff[3] ** 2 * es[3][3])

        if which == 'st_coupling':
            energy = 0.5 * (diff[2] * diff[3] * es[2][3])

        if which == 'bs_coupling':
            energy = 0.5 * ((diff[0] * diff[2] * es[0][2])
                            + (diff[1] * diff[2] * es[1][2]))

        if which == 'bt_coupling':
            energy = 0.5 * ((diff[0] * diff[3] * es[0][3])
                            + (diff[1] * diff[3] * es[1][3]))

        if which == 'bb_coupling':
            energy = 0.5 * (diff[0] * diff[1] * es[0][1])

        if which == 'st':
            energy = 0.5 * ((diff[0] ** 2 * es[0][0])
                            + (diff[1] ** 2 * es[1][1])
                            + (diff[2] ** 2 * es[2][2])
                            + (diff[3] ** 2 * es[3][3])
                            + (diff[2] * diff[3] * es[2][3]))

        if which == 'bs':
            energy = 0.5 * ((diff[0] ** 2 * es[0][0])
                            + (diff[1] ** 2 * es[1][1])
                            + (diff[2] ** 2 * es[2][2])
                            + (diff[0] * diff[2] * es[0][2])
                            + (diff[1] * diff[2] * es[1][2]))

        if which == 'bt':
            energy = 0.5 * ((diff[0] ** 2 * es[0][0])
                            + (diff[1] ** 2 * es[1][1])
                            + (diff[3] ** 2 * es[3][3])
                            + (diff[0] * diff[3] * es[0][3])
                            + (diff[1] * diff[3] * es[1][3]))

        return energy

    def calculateLocalElasticity(self, bp, frames=None, helical=False, unit='kT'):
        r"""Calculate local elastic matrix or stiffness matrix

        In case of :ref:`base-step-image`: Shift (:math:`Dx`), Slide (:math:`Dy`), Rise (:math:`Dz`),
        Tilt (:math:`\tau`), Roll (:math:`\rho`) and Twist (:math:`\omega`), following elastic matrix is calculated.

        .. math::

            \mathbf{K}_{base-step} = \begin{bmatrix}
            K_{Dx}        & K_{Dx,Dy}      & K_{Dx,Dz}      & K_{Dx,\tau}      & K_{Dx,\rho}      & K_{Dx,\omega} \\
            K_{Dx,Dy}     & K_{Dy}         & K_{Dy,Dz}      & K_{Dy,\tau}      & K_{Dy,\rho}      & K_{Dy,\omega} \\
            K_{Dx,Dz}     & K_{Dy,Dz}      & K_{Dz}         & K_{Dz,\tau}      & K_{Dz,\rho}      & K_{Dz,\omega} \\
            K_{Dx,\tau}   & K_{Dy,\tau}    & K_{Dz,\tau}    & K_{\tau}         & K_{\tau, \rho}   & K_{\tau,\omega} \\
            K_{Dx,\rho}   & K_{Dy,\rho}    & K_{Dz,\rho}    & K_{\tau, \rho}   & K_{\rho}         & K_{\rho,\omega} \\
            K_{Dx,\omega} & K_{Dy,\omega}  & K_{Dz,\omega}  & K_{\tau, \omega} & K_{\rho, \omega} & K_{\omega} \\
            \end{bmatrix}


        In case of :ref:`helical-base-step-image`: x-displacement (:math:`dx`), y-displacement (:math:`dy`), h-rise (:math:`h`),
        inclination (:math:`\eta`), tip (:math:`\theta`) and twist (:math:`\Omega`), following elastic matrix is calculated.

        .. math::

            \mathbf{K}_{helical-base-step} = \begin{bmatrix}
            K_{dx}        & K_{dx,dy}      & K_{dx,h}      & K_{dx,\eta}      & K_{dx,\theta}      & K_{dx,\Omega} \\
            K_{dx,dy}     & K_{dy}         & K_{dy,h}      & K_{dy,\eta}      & K_{dy,\theta}      & K_{dy,\Omega} \\
            K_{dx,h}      & K_{dy,h}       & K_{h}         & K_{h,\eta}       & K_{h,\theta}       & K_{h,\Omega} \\
            K_{dx,\eta}   & K_{dy,\eta}    & K_{h,\eta}    & K_{\eta}         & K_{\eta, \theta}   & K_{\eta,\Omega} \\
            K_{dx,\theta} & K_{dy,\theta}  & K_{h,\theta}  & K_{\eta, \theta} & K_{\theta}         & K_{\theta,\Omega} \\
            K_{dx,\Omega} & K_{dy,\Omega}  & K_{h,\Omega}  & K_{\eta, \Omega} & K_{\theta, \Omega} & K_{\Omega} \\
            \end{bmatrix}


        Parameters
        ----------
        bp : list
            List of two base-steps forming the DNA segment.
            For example: with ``bp=[5, 50]``, 5-50 base-step segment will be considered.

        frames : list
            List of two trajectory frames between which parameters will be extracted. It can be used to select portions
            of the trajectory. For example, with ``frames=[100, 1000]``, 100th to 1000th frame of the trajectory will be
            considered.

        helical : bool
            If ``helical=True``, elastic matrix for **helical base-step** parameters are calculated. Otherwise,
            by default, elastic matrix for **base-step** parameters are calculated.

        unit : str
            Unit of energy. Allowed units are: ``'kT', 'kJ/mol' and 'kcal/mol'``.

        Return
        ------
        mean : numpy.ndarray
            Value of parameters at which energy is zero. Minimum point on energy landscape.

            if ``helical=False``

            .. math::
                \begin{bmatrix}
                    Dx_0  &  Dy_0 & Dz_0 & \tau_0 & \rho_0 & \omega_0
                \end{bmatrix}

            if ``helical=True``

             .. math::
                \begin{bmatrix}
                    dx_0  &  dy_0 & h_0 & \eta_0 & \theta_0 & \Omega_0
                \end{bmatrix}

        result : numpy.ndarray
            Elastic matrix.

        """
        acceptedUnit = ['kT', 'kJ/mol', 'kcal/mol']

        if unit not in acceptedUnit:
            raise ValueError(" {0} not accepted. Use any of the following: {1} ".format(unit, acceptedUnit))

        frames = self._validateFrames(frames)

        name = '{0}-{1}-{2}-{3}-local-{4}'.format(bp[0], bp[1], frames[0], frames[1], int(helical))


        if name not in self.esMatrix:
            time, array = self._extractLocalParameters(self.dna, bp, helical=helical, frames=frames)
            mean = np.mean(array, axis = 1)
            esMatrix = self.getElasticMatrix(array)
            self.esMatrix[name] = esMatrix
            self.minimumPoint[name] = mean
        else:
            esMatrix = self.esMatrix[name]
            mean = self.minimumPoint[name]

        if unit == 'kJ/mol':
            result = 2.4946938107879997 * esMatrix      # (1.38064852e-23 * 300 * 6.023e23 / 1000 )  kT.NA/1000
        elif unit == 'kcal/mol':
            result = 0.5962461306854684 * esMatrix      # (1.38064852e-23 * 300 * 6.023e23 / 1000 / 4.184)  kT.NA/1000
        else:
            result = esMatrix

        return mean, result

    def getLocalElasticityByTime(self, bp, skip, helical=False, unit='kT'):
        r"""Calculate local elastic properties as a function of time for convergence check

        It can be used to obtained elastic properties as a function of time.

        .. note:: Elastic properties cannot be calculated using a single frame because fluctuations are required.
                  Therefore, here time means trajectory between zero time to given time.


        When ``helical='False'``, following is obtained:

            1) Shift (:math:`K_{Dx}`)
            2) Slide (:math:`K_{Dy}`)
            3) Rise (:math:`K_{Dz}`)
            4) Tilt (:math:`K_{\tau}`)
            5) Roll (:math:`K_{\rho}`)
            6) Twist (:math:`K_{\omega}`)
            7) :math:`K_{Dx,Dy}`
            8) :math:`K_{Dy,Dz}`
            9) :math:`K_{Dz,\tau}`
            10) :math:`K_{\tau, \rho}`
            11) :math:`K_{\rho,\omega}`
            12) :math:`K_{Dx,Dz}`
            13) :math:`K_{Dy,\tau}`
            14) :math:`K_{Dz,\rho}`
            15) :math:`K_{\tau,\omega}`
            16) :math:`K_{Dx,\tau}`
            17) :math:`K_{Dy,\rho}`
            18) :math:`K_{Dz,\omega}`
            19) :math:`K_{Dx,\rho}`
            20) :math:`K_{Dy,\omega}`
            21) :math:`K_{Dx,\omega}`

        When ``helical='True'``, following is obtained:

            1) Shift (:math:`K_{Dx}`)
            2) Slide (:math:`K_{Dy}`)
            3) Rise (:math:`K_{h}`)
            4) Tilt (:math:`K_{\eta}`)
            5) Roll (:math:`K_{\theta}`)
            6) Twist (:math:`K_{\Omega}`)
            7) :math:`K_{dx,dy}`
            8) :math:`K_{dy,h}`
            9) :math:`K_{h,\eta}`
            10) :math:`K_{\eta, \theta}`
            11) :math:`K_{\theta,\Omega}`
            12) :math:`K_{dx,h}`
            13) :math:`K_{dy,\eta}`
            14) :math:`K_{h,\theta}`
            15) :math:`K_{\tau,\Omega}`
            16) :math:`K_{dx,\eta}`
            17) :math:`K_{dy,\theta}`
            18) :math:`K_{h,\Omega}`
            19) :math:`K_{dx,\theta}`
            20) :math:`K_{dy,\Omega}`
            21) :math:`K_{dx,\Omega}`

        .. currentmodule:: dnaMD

        Parameters
        ----------
        bp : list
            List of two base-steps forming the DNA segment.
            For example: with ``bp=[5, 50]``, 5-50 base-step segment will be considered.

        skip : int
            How many frames to skip for next time-frame. Lower the number, slower will be the calculation.

        helical : bool
            If ``helical=True``, elastic matrix for **helical base-step** parameters are calculated. Otherwise,
            by default, elastic matrix for **base-step** parameters are calculated.

        unit : str
            Unit of energy. Allowed units are: ``'kT', 'kJ/mol' and 'kcal/mol'``.


        Returns
        -------
        time : numpy.ndarray
            1D array containing time values.

        modulus : numpy.ndarray
            2D array of shape (21, nframes). Order of elasticity in this array is same as to the above listed order.

        """

        length = len(self.dna.time[:])

        time, elasticity = [], []
        for i in range(skip, length, skip):

            props = None

            mean, esy_t = self.calculateLocalElasticity(bp, helical=helical, unit=unit)

            props = np.diagonal(esy_t, offset=0)
            props = np.hstack((props, np.diagonal(esy_t, offset=1)))
            props = np.hstack((props, np.diagonal(esy_t, offset=2)))
            props = np.hstack((props, np.diagonal(esy_t, offset=3)))
            props = np.hstack((props, np.diagonal(esy_t, offset=4)))
            props = np.hstack((props, np.diagonal(esy_t, offset=5)))

            time.append(self.dna.time[i])
            elasticity.append(props)

        elasticity = np.asarray(elasticity)

        return time, elasticity.T

    def getLocalDeformationEnergy(self, bp, complexDna, freeDnaFrames=None, boundDnaFrames=None, helical=False, unit='kT', which='all'):
        r"""Deformation energy of the input DNA using local elastic properties

        It can be used to calculate deformation energy of a base-step/s from input DNA object with reference to
        the same base-step/s DNA present in the current DNA object.

        The deformation free energy is calculated using elastic matrix as follows

        .. math::

            G = \frac{1}{2}\mathbf{xKx^T}

        When ``helical='False'``

        .. math::
            \mathbf{K} = \mathbf{K}_{base-step}

        .. math::

            \mathbf{x} =  \begin{bmatrix}
                              (Dx_{i}-Dx_0)  &  (Dy_i - Dy_0) & (Dz_i - Dz_0) & (\tau_i - \tau_0) &
                              (\rho_i - \rho_0) & (\omega_i - \omega_0)
                          \end{bmatrix}


        When ``helical='True'``

        .. math::
            \mathbf{K} = \mathbf{K}_{helical-base-step}

        .. math::
            \mathbf{x} =  \begin{bmatrix}
                              (dx_{i}-dx_0)  &  (dy_i - dy_0) & (h_i - h_0) & (\eta_i - \eta_0) &
                              (\theta_i - \theta_0) & (\Omega_i - \Omega_0)
                          \end{bmatrix}

        .. currentmodule:: dnaMD


        Parameters
        ----------
        bp : list
            List of two base-steps forming the DNA segment.
            For example: with ``bp=[5, 50]``, 5-50 base-step segment will be considered.
        complexDna : :class:`dnaMD.DNA`
            Input :class:`dnaMD.DNA` instance for which deformation energy will be calculated.
        freeDnaFrames : list
            To select a trajectory segment of current (free) DNA data.
            List of two trajectory frames between which parameters will be extracted. It can be used to select portions
            of the trajectory. For example, with ``frames=[100, 1000]``, 100th to 1000th frame of the trajectory will be
            considered.
        boundDnaFrames : list
            To select a trajectory segment of input (bound) DNA data.
            List of two trajectory frames between which parameters will be extracted. It can be used to select portions
            of the trajectory. For example, with ``frames=[100, 1000]``, 100th to 1000th frame of the trajectory will be
            considered.
        helical : bool
            If ``helical=True``, elastic matrix for **helical base-step** parameters are calculated. Otherwise,
            by default, elastic matrix for **base-step** parameters are calculated.
        unit : str
            Unit of energy. Allowed units are: ``'kT', 'kJ/mol' and 'kcal/mol'``.
        which : all
            For which motions, energy should be calculated.

            Following keywords are available:
                * ``'all'`` : Use entire elastic matrix -- all parameters with their coupling
                * ``'diag'`` : Use diagonal of elastic matrix -- all motions but no coupling
                * ``'shift'`` or ``'x-disp'``
                * ``'slide'`` or ``'y-idsp'``
                * ``'rise'`` or ``'h-rise'``
                * ``'tilt'`` or ``'inclination'``
                * ``'roll'`` or ``'tip'``
                * ``'twist'`` or ``'h-twist'``

        Returns
        -------
        time : numpy.ndarray
            1D array containing time values.

        energy : numpy.ndarray
            1D array of shape (nframes) containing energy value for query DNA.

        """
        if which not in self.enLocalTypes:
            raise ValueError('{0} is not a supported energy keywords.\n Use any of the following: \n {1}'.format(
                which, self.enLocalTypes))

        if helical and which not in ['all', 'diag', 'x-disp', 'y-disp', 'h-rise', 'inclination', 'tip', 'h-twist']:
            raise ValueError (' {0} is not a helical base-step parameters'.format(which))

        if not helical and which not in ['all', 'diag', 'shift', 'slide', 'rise', 'tilt', 'roll', 'twist']:
            raise ValueError (' {0} is not a base-step parameters'.format(which))

        means, esMatrix = self.calculateLocalElasticity(bp, frames=freeDnaFrames, helical=helical, unit=unit)

        time, array = self._extractLocalParameters(complexDna, bp, frames=boundDnaFrames, helical=helical)

        energy = []
        diffs = []
        for i in range(array[0].shape[0]):
            vec = array[:, i]
            diff = vec - means
            diffs.append(diff)
            t_energy = self._calcLocalEnergy(diff, esMatrix, which)
            energy.append(t_energy)

        return time, np.asarray(energy)

    def _calcLocalEnergy(self, diff, es, which):
        r"""Calculate local deformation energy using a difference vector.

        It is called in :meth:`dnaEY.getLocalDeformationEnergy` for energy calculation of each frame.

        Parameters
        ----------
        diff : numpy.ndarray
            Array of difference between minimum and current parameter values.

        es : numpy.ndarray
            Elastic matrix. See in :meth:`dnaEY.calculateLocalElasticity` about elastic matrix.

        which : str
            For which type of motions, energy will be calculated.
            see ``which`` parameter in :meth:`dnaEY.getLocalDeformationEnergy` for keywords.

        Return
        ------
        energy : float
            Deformation free energy value

        """

        if which not in self.enLocalTypes:
            raise ValueError('{0} is not a supported energy keywords.\n Use any of the following: \n {1}'.format(
                which, self.enLocalTypess))

        energy = None

        if which == 'all':
            temp = np.matrix(diff)
            energy = 0.5 * ((temp * es) * temp.T)
            energy = energy[0,0]

        if which == 'diag':
            energy = 0.5 * ((diff[0] ** 2 * es[0][0])
                            + (diff[1] ** 2 * es[1][1])
                            + (diff[2] ** 2 * es[2][2])
                            + (diff[3] ** 2 * es[3][3])
                            + (diff[4] ** 2 * es[4][4])
                            + (diff[5] ** 2 * es[5][5]))

        if which == 'shift' or which == 'x-disp':
            energy = 0.5 * (diff[0] ** 2 * es[0][0])

        if which == 'slide' or which == 'y-disp':
            energy = 0.5 * (diff[1] ** 2 * es[1][1])

        if which == 'rise' or which == 'h-rise':
            energy = 0.5 * (diff[2] ** 2 * es[2][2])

        if which == 'tilt' or which == 'inclination':
            energy = 0.5 * (diff[3] ** 2 * es[3][3])

        if which == 'roll' or which == 'tip':
            energy = 0.5 * (diff[4] ** 2 * es[4][4])

        if which == 'twist' or which == 'h-twist':
            energy = 0.5 * (diff[5] ** 2 * es[5][5])

        return energy
