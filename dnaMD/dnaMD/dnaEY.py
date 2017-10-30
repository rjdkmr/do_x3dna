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

from . import dnaMD


class eyDNA:

    def __init__(self, num_bp, filename=None, startBP=1):
        self.dna = dnaMD.DNA(num_bp, filename=filename, startBP=startBP)
        self.esMatrix = dict()
        self.minimumPoint = dict()

    def getElasticMatrix(self, inputArray):
        inputArray = np.array( inputArray )

        # Calculation of covariance matrix
        CovMat = np.cov(inputArray, bias=1)

        # Change to a matrix object
        CovMat = np.matrix(CovMat)

        # Inverse of the covariance matrix
        InvCovMat = CovMat.I

        return np.asarray(InvCovMat)


    def getStretchTwistBend(self, bp, frames=None, paxis='Z', masked=True):

        if frames is None:
            frames = [0, -1]
        else:
            if (len(frames) != 2):
                raise ValueError("frames should be a list containing lower and higher limit. See, documentation!!!")

            if frames[1] != -1 and frames[0] > frames[1]:
                raise ValueError("frames should be a list containing lower and higher limit. See, documentation!!!")

        if (len(bp) != 2):
            raise ValueError("bp should be a list containing first and last bp of a segment. See, documentation!!!")

        if bp[0] > bp[1]:
            raise ValueError("bp should be a list containing first and last bp of a segment. See, documentation!!!")


        time, clen = self.dna.time_vs_parameter('h-rise', bp=bp, merge=True, merge_method='sum', masked=masked)
        clen = np.asarray(clen) * 0.1  # conversion to nm

        time, htwist = self.dna.time_vs_parameter('h-twist', bp=bp, merge=True, merge_method='sum', masked=masked)
        htwist = np.deg2rad(htwist)  # Conversion to radian

        angleOne, angleTwo = self.dna.calculate_2D_angles_bw_tangents(paxis, bp, masked=masked)


        # Rarely there are nan during angle calculation, remove those nan
        nanInOne = np.isnan( angleOne[frames[0]:frames[1]] )
        nanInTwo = np.isnan( angleTwo[frames[0]:frames[1]] )
        notNan = ~(nanInOne + nanInTwo)
        notNanIdx = np.nonzero( notNan )

        if frames[1] == -1:
            array = np.array( [ angleOne[frames[0]:][notNanIdx], angleTwo[frames[0]:][notNanIdx],
                                clen[frames[0]:][notNanIdx], htwist[frames[0]:][notNanIdx] ] )
        else:

            array = np.array( [ angleOne[frames[0]:frames[1]][notNanIdx], angleTwo[frames[0]:frames[1]][notNanIdx],
                                clen[frames[0]:frames[1]][notNanIdx], htwist[frames[0]:frames[1]][notNanIdx] ] )

        mean = np.mean(array, axis = 1)

        esMatrix = self.getElasticMatrix(array)

        modulus = 4.1419464 * np.array(esMatrix) * mean[2]

        return mean, modulus

    def getStretchTwist(self, bp, frames=None, masked=False):

        if frames is None:
            frames = [0, -1]
        else:
            if (len(frames) != 2):
                raise ValueError("frames should be a list containing lower and higher limit. See, documentation!!!")

            if frames[1] != -1 and frames[0] > frames[1]:
                raise ValueError("frames should be a list containing lower and higher limit. See, documentation!!!")

        if (len(bp) != 2):
            raise ValueError("bp should be a list containing first and last bp of a segment. See, documentation!!!")

        if bp[0] > bp[1]:
            raise ValueError("bp should be a list containing first and last bp of a segment. See, documentation!!!")

        time, clen = self.dna.time_vs_parameter('h-rise', bp=bp, merge=True, merge_method='sum', masked=masked)
        clen = np.asarray(clen) * 0.1  # conversion to nm

        time, htwist = self.dna.time_vs_parameter('h-twist', bp=bp, merge=True, merge_method='sum', masked=masked)
        htwist = np.deg2rad(htwist)  # Conversion to radian

        if frames[1] == -1:
            array = np.array( [clen[frames[0]:], htwist[frames[0]:] ] )
        else:
            array = np.array( [ clen[frames[0]:frames[1]], htwist[frames[0]:frames[1]]] )

        mean = np.mean(array, axis = 1)

        esMatrix = self.getElasticMatrix(array)

        modulus = 4.1419464 * np.array(esMatrix) * mean[0]

        return mean, modulus

    def getElasticityByTime(self, esType, bp, skip, masked=False, paxis='Z'):
        """

            1) bending-1

            2) bending-2

            3) Stretching

            4) Twisting

            5) bending-1-bending-2

            6) bending-2-stretching

            7) Stretching-Twisting

            8) bending-1-stretching

            9) bending2-Twisting

            10) bending-1-twisting

        """

        if esType not in ['BST', 'ST']:
            raise ValueError('Accepted keywords are BST and ST !!!')

        length = len(self.dna.time[:])

        time, modulus = [], []
        for i in range(skip, length, skip):

            props = None
            if esType == 'BST':
                mean, modulus_t = self.getStretchTwistBend(bp, frames=[0, i], paxis=paxis, masked=True)

                props = np.diagonal(modulus_t, offset=0)
                props = np.hstack((props, np.diagonal(modulus_t, offset=1)))
                props = np.hstack((props, np.diagonal(modulus_t, offset=2)))
                props = np.hstack((props, np.diagonal(modulus_t, offset=3)))

            if esType == 'ST':
                mean, modulus_t = self.getStretchTwist(bp, frames=[0, i], masked=masked)

                props = np.diagonal(modulus_t, offset=0)
                props = np.hstack((props, np.diagonal(modulus_t, offset=1)))

            time.append(self.dna.time[i])
            modulus.append(props)

        modulus = np.asarray(modulus)

        return time, modulus.T


def calc_energy_by_axis(RefDna, SubjDna, bp, axis ='Z', bp_range=True, windows=10, dof=['h-Twist','h-Rise'], err_type='acf'):

    RefAxis, Ref_bp_idx = SubjDna.get_parameters('Helical {0}-axis' .format(axis),bp, bp_range=bp_range)

    mean_axis = np.mean(RefAxis,axis=1)
    maxAxis = np.amax(mean_axis)
    minAxis = np.amin(mean_axis)
    axis_range = (maxAxis-minAxis)/windows

    bp_index = []
    final_axis = []
    for i in range(windows):
        start = minAxis+(i*axis_range)
        end = start + axis_range
        idx = []
        for j in range(len(mean_axis)):
            if((start <= mean_axis[j]) and (end>mean_axis[j])):
                idx.append(Ref_bp_idx[j])
        if(len(idx)>0):
            final_axis.append( start + (end-start)/2 )
            bp_index.append(idx)

    energy = []
    error = []
    for i in range(len(final_axis)):
        tmp_en, tmp_err = calc_deform_energy(RefDna, SubjDna, bp_index[i], bp_range=False, dof=dof, err_type=err_type)
        energy.append(tmp_en)
        error.append(tmp_err)

    return final_axis, energy, error

def calc_energy_by_base_steps(RefDna, SubjDna, bp, bp_range=True, dof=['h-Twist','h-Rise'], merge_bp=1, err_type='acf'):

    if(not bp_range) and (merge_bp>1):
        print ("\nERROR: Merging of base pairs/steps only possible with given base pair/steps range\n")
        exit(1)

    length = bp[1]-bp[0]+1
    bins = int(length/merge_bp)

    mid_bin = []
    energy = []
    error = []

    i = bp[0]
    while(1):
        start = i
        if((i + merge_bp ) >= length):
            end = i + (length-i)
        else:
            end   = i + merge_bp

        if(start < end):
            mid_bin.append(int(start+(end-start)/2))
            tmp_en, tmp_err = calc_deform_energy(RefDna, SubjDna, [start, end], bp_range=True, dof=dof, err_type=err_type)
            energy.append(tmp_en)
            error.append(tmp_err)

        if(i >= length):
            break

        i += merge_bp

    return mid_bin, energy, error


def calc_energy_landscape(dna, bp, range_limit, dof=['h-Twist','h-Rise'], window=[50,50], bp_range=True):

    MeanParameters, elastic_constant = calc_force_constant(dna, bp, bp_range=bp_range, dof=dof)

    binx = range_limit[0]/window[0]
    biny = range_limit[1]/window[1]

    DataX = np.arange(MeanParameters[0]-range_limit[0], MeanParameters[0]+range_limit[0], binx)
    DataY = np.arange(MeanParameters[1]-range_limit[1], MeanParameters[1]+range_limit[1], biny)

    deformation_x = np.subtract(DataX, MeanParameters[0])
    deformation_y = np.subtract(DataY, MeanParameters[1])

    deformation = []
    energy = []
    for i in range(len(deformation_x)):
        tmp_energy = []
        for j in range(len(deformation_y)):
            mat = np.array([deformation_x[i], deformation_y[j]])
            tmpen = 0.5 * np.dot(mat, np.dot(elastic_constant,mat.T))
            tmp_energy.append(tmpen)
        energy.append(tmp_energy)

    energy = np.array(energy)

    return DataX, DataY, energy

def calc_deform_energy(RefDna, SubjDna, bp, bp_range=True, dof=['h-Twist','h-Rise'], err_type='acf'):

    MeanRefParameters, elastic_constant = calc_force_constant(RefDna, bp, bp_range=bp_range, dof=dof)

    SubjParameters = []
    for i in range(len(dof)):
        parameter, bp_idx = SubjDna.get_parameters(dof[i], bp, bp_range=bp_range)
        SubjParameters.append(parameter)

    Sum_SubjParameters = []
    for i in range(len(dof)):
        Sum_SubjParameters.append( np.sum(SubjParameters[i], axis=0) )

    dev_parameters = []
    for i in range(len(dof)):
        dev_parameters.append(np.mean(Sum_SubjParameters[i]) - MeanRefParameters[i])

    mat = np.array(dev_parameters)
    energy = 0.5 * np.dot(mat, np.dot(elastic_constant,mat.T))

    all_energy = []
    for i in range(len(Sum_SubjParameters[0])):
        dev_parameters = []
        for j in range(len(dof)):
            dev_parameters.append(Sum_SubjParameters[j][i] - MeanRefParameters[j])

        mat = np.array(dev_parameters)
        tmp_energy = 0.5 * np.dot(mat, np.dot(elastic_constant,mat.T))
        all_energy.append(tmp_energy)
    energy_error = (dm.get_error(SubjDna.time,[all_energy],1,err_type=err_type))[0]

    return energy, energy_error

def calc_elastic_constants(dna, bp, bp_range=True, dof=['h-Twist','h-Rise']):

    parameters = []
    for i in range(len(dof)):
        parameter, bp_idx = dna.get_parameters(dof[i], bp, bp_range=bp_range)
        parameters.append(parameter)

    sum_parameters = []
    for i in range(len(parameters)):
        sum_parameters.append( np.sum(parameters[i], axis=0) )

    Mean_parameters = []
    for i in range(len(parameters)):
        Mean_parameters.append( np.mean(sum_parameters[i]) )

    print("Mean values =\n", Mean_parameters)

    # Calculation of covariance matrix
    array = np.array(sum_parameters)
    CovMat = np.cov(array,bias=1)

    # Change to a matrix object
    CovMat = np.matrix(CovMat)

    # Inverse of the covariance matrix
    InvCovMat = CovMat.I

    print("Spring constant (kT/A or rad) =\n", InvCovMat )

    # Conversion of matrix to array
    InvCovMat = np.array(InvCovMat)

    return Mean_parameters, InvCovMat
