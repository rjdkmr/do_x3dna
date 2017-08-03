#!/usr/bin/env python

#AUTHOR: Rajendra Kumar
#EMAIL: rkumar@gwdg.de
#DATE: Jun 27, 2013
#ARCH: any unix
#REQUIRES: numpy, dna_main
#STATUSCOMMENTS: Version 1.0
#DESCRIPTION: To calculate elastic force constant and deformation energy of the DNA
#USAGE: import dna_elasticity
#EXAMPLE:
###ENDHEADER###


import numpy as np
import dnaMD as dm

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

	# Convertion of matrix to array
	InvCovMat = np.array(InvCovMat)

	return Mean_parameters, InvCovMat

def calc_ST_props(dna, bp, dof=['h-Rise','h-Twist'], frames=[0, -1], masked=False):

	allowed_parms = ['h-Rise','h-Twist', 'Rise', 'Twist']
	unit_conv = { 'h-Rise':0.1, 'h-Twist': (np.pi/180), 'Rise':0.1, 'Twist': (np.pi/180)}

	for p in range(len(dof)):
		if not dof[p] in allowed_parms:
			raise ValueError('\nParameter {0} is not allowed... Allowed parameters: {1}' .format(dof[p], allowed_parms))

	if (len(frames) != 2):
		raise ValueError("See, documentation for frames usage!!!")

	if frames[1] != -1 and frames[0] > frames[1]:
		raise ValueError("See, documentation for frames usage!!!")

	if (len(bp) != 2):
		raise ValueError("See, documentation for bp usage!!!")

	if bp[0] > bp[1]:
		raise ValueError("See, documentation for bp usage!!!")

	parameters = []
	rise_idx = -9999
	for i in range(len(dof)):
		time, parameter = dna.time_vs_parameter(dof[i], bp, merge=True, merge_method='sum', masked=masked)

		if dof[i]== 'h-Rise' or dof[i]== 'Rise':
			rise_idx = i

		if frames[1] == -1:
			parameters.append(np.multiply(parameter[frames[0]:], unit_conv[dof[i]]))
		else:
			parameters.append(np.multiply(parameter[frames[0]:frames[1]], unit_conv[dof[i]]))

	mean = []
	for parameter in parameters:
		mean.append(np.mean(parameter))

	# Calculation of covariance matrix
	array = np.array(parameters)
	CovMat = np.cov(array,bias=1)

	# Change to a matrix object
	CovMat = np.matrix(CovMat)

	# Inverse of the covariance matrix
	InvCovMat = CovMat.I

	# Convertion of matrix to array
	InvCovMat = np.array(InvCovMat)

	modulus = 4.1419464 * mean[rise_idx] * InvCovMat

	return mean, CovMat, modulus


def calc_STB_modulus(dna, bp, frames=[0, -1], paxis='Z', masked=True):


	if (len(frames) != 2):
		raise ValueError("See, documentation for frames usage!!!")

	if frames[1] != -1 and frames[0] > frames[1]:
		raise ValueError("See, documentation for frames usage!!!")

	if (len(bp) != 2):
		raise ValueError("See, documentation for bp usage!!!")

	if bp[0] > bp[1]:
		raise ValueError("See, documentation for bp usage!!!")


	time, clen = dna.time_vs_parameter('h-Rise', bp=bp, merge=True, merge_method='sum', masked=True)
	clen = np.asarray(clen) * 0.1  # conversion to nm

	time, htwist = dna.time_vs_parameter('h-Twist', bp=bp, merge=True, merge_method='sum', masked=True)
	htwist = np.deg2rad(htwist)  # Conversion to radian

	angleOne, angleTwo = dna.calculate_2D_angles_bw_tangents(paxis, bp)


	if frames[1] == -1:
		# Rarely there are nan during angle calculation, remove those nan
		nanInOne = np.isnan( angleOne[frames[0]:] )
		nanInTwo = np.isnan( angleTwo[frames[0]:] )
		notNan = ~(nanInOne + nanInTwo)
		notNanIdx = np.nonzero( notNan )

		array = np.array( [ angleOne[frames[0]:][notNanIdx], angleTwo[frames[0]:][notNanIdx],
							clen[frames[0]:][notNanIdx], htwist[frames[0]:][notNanIdx] ] )
	else:
		# Rarely there are nan during angle calculation, remove those nan
		nanInOne = np.isnan( angleOne[frames[0]:frames[1]] )
		nanInTwo = np.isnan( angleTwo[frames[0]:frames[1]] )
		notNan = ~(nanInOne + nanInTwo)
		notNanIdx = np.nonzero( notNan )

		array = np.array( [ angleOne[frames[0]:frames[1]][notNanIdx], angleTwo[frames[0]:frames[1]][notNanIdx],
							clen[frames[0]:frames[1]][notNanIdx], htwist[frames[0]:frames[1]][notNanIdx] ] )

	mean = np.mean(array, axis = 1)

	# Calculation of covariance matrix
	CovMat = np.cov(array,bias=1)

	# Change to a matrix object
	CovMat = np.matrix(CovMat)

	# Inverse of the covariance matrix
	InvCovMat = CovMat.I

	modulus = 4.1419464 * np.array(InvCovMat) * mean[2]

	return mean, InvCovMat, modulus
