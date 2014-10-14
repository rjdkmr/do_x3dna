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
####### VIM Setting ###########
#
# set smarttab
# set tabstop=4
# set shiftwidth=4
# set autoindent
# set number
#
###############################

import numpy as np
import re, os, sys
import string, random
import math
import subprocess as sub
try:
	from scipy.interpolate import splprep, splev
	scipy_imported = True
except:
	scipy_imported = False

class base_step:
	def __init__(self):
		self.shift = []
		self.slide = []
		self.rise = []
		self.tilt = []
		self.roll = []
		self.twist = []

		self.x_disp = []
		self.y_disp = []
		self.h_rise = []
		self.inclination = []
		self.tip = []
		self.h_twist = []

		self.major_pp = []
		self.minor_pp = []
		self.major_refine = []
		self.minor_refine = []

		self.hel_Xaxis = []
		self.hel_Yaxis = []
		self.hel_Zaxis = []
		
		self.hel_Xaxis_smth = []
		self.hel_Yaxis_smth = []
		self.hel_Zaxis_smth = []

		self.curvature = []
		self.tangent = []

class base_pair():
	def __init__(self):
		self.shear = []
		self.stretch = []
		self.stagger = []
		self.buckle = []
		self.propeller = []
		self.opening = []
		
		self.hbond = []

		self.radS1 = []
		self.radS2 = []

		self.alpha_s1 = []
		self.beta_s1 = []
		self.gamma_s1 = []
		self.delta_s1 = []
		self.epsilon_s1 = []
		self.zeta_s1 = []
		self.chi_s1 = []
		
		self.alpha_s2 = []
		self.beta_s2 = []
		self.gamma_s2 = []
		self.delta_s2 = []
		self.epsilon_s2 = []
		self.zeta_s2 = []
		self.chi_s2 = []

class DNA:
	"""DNA class stores all data obtained from the input files.
		
		**To initialize this class:** ::
		
				dna = DNA(60)       # 60 is the number of basepairs

		This class also contains several methods (functions) that are discussed in following sections.

	"""

	def __init__(self,num_bp):
		self.num_bp = num_bp
		self.num_step = num_bp-1
		self.smooth_axis = False
		self.time = []
		self.mask = []
		self.base_steps = []
		self.base_pairs = []

		for i in range(num_bp):
			self.base_pairs.append(base_pair())

		for i in range(self.num_step):
			self.base_steps.append(base_step())

	def get_parameters(self,parameter, bp, bp_range=True, masked=False):
		"""To get the parameters over all frame for the given range of base pair/steps
		
		**Arguments:**
			- ``parameter (string)``: [Name of the prameter]
				Currently accepted keywords are as follows:
					* ``Shear``
					* ``Stretch``
					* ``Stagger``
					* ``Buckle``
					* ``Propeller``
					* ``Opening``
					* ``Shift``
					* ``Slide``
					* ``Rise``
					* ``Tilt``
					* ``Roll``
					* ``Twist``
					* ``X-disp``
					* ``Y-disp``
					* ``h-Rise``
					* ``Inclination``
					* ``Tip``
					* ``h-Twist``
					* ``Helical X-axis``
					* ``Helical Y-axis``
					* ``Helical Z-axis``
					* ``Helical X-axis smooth``
					* ``Helical Y-axis smooth``
					* ``Helical Z-axis smooth``
					* ``Helical axis curvature``
					* ``Helical axis tangent``
					* ``Radius S-1``
					* ``Radius S-2``
					* ``Major Groove``
					* ``Major Groove Refined``
					* ``Minor Groove``
					* ``Minor Groove Refined``
					* ``alpha S-1``
					* ``beta S-1``
					* ``gamma S-1``
					* ``delta S-1``
					* ``epsilon S-1``
					* ``zeta S-1``
					* ``chi S-1``
					* ``alpha S-2``
					* ``beta S-2``
					* ``gamma S-2``
					* ``delta S-2``
					* ``epsilon S-2``
					* ``zeta S-2``
					* ``chi S-2``

			- ``bp (1D list) or (1D array)``: base-pairs to analyze
				Example: ::

							bp = [6]                                # bp_range = False
							bp = [4,15]                             # bp_range = True
							bp = range(4,15)                        # bp_range = False
							bp = np.arange(4,15)                    # bp_range = False
							bp = [2,5,6,7,9,12,18]                  # bp_range = False
					
			* ``bp_range (bool)``: ``Dfault=True``: As shown above, if ``True``, bp is taken as a range otherwise list or numpy array
			
			
			* ``masked (bool)``: ``Dfault=False``: To skip specific frames/snapshots. dnaMD.DNA.mask array should be set to use this functionality. This array contains boolean (either ``True`` or ``False``) value for each frame to mask the frames. Presently, mask array is automatically generated during :meth:`dnaMD.DNA.generate_smooth_axis` method to skip those frames where 3D fitting curve was not successfull within the given critera.
            

        **Returns:**
            - ``parameters[bp][nframe] (2D list)``: where bp is number of base pairs/steps and nframe is total number of frames in the trajectory.
		
		"""

		bp_idx, dum = get_idx_of_bp_parameters(bp,[],bp_range)
		append = False
		empty = False
		key = 'dummy'
		idx = 0
		data=[]
		midx = []


		# Masking values according to mask array
		if masked and len(self.mask)==0:
			raise ValueError("mask array is not set. mask array is set within generate_smooth_axis() \n")
		
		for i in range(len(self.time)):
			if masked:
				if	self.mask[i] == False:
					midx.append(i)
			else:
				midx.append(i)

		#Extracting data for given base pairs and parameters combination
		#Base pair parameters
		if(parameter=='Shear'):
			key = 'Shear'
			for i in range(len(bp_idx)):
				if (len(self.base_pairs[bp_idx[i]].shear) == 0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_pairs[bp_idx[i]].shear[midx])
					append =True

		if(parameter=='Stretch'):
			key = 'Stretch'
			for i in range(len(bp_idx)):
				if (len(self.base_pairs[bp_idx[i]].stretch)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_pairs[bp_idx[i]].stretch[midx])
					append = True

		if(parameter=='Stagger'):
			key = 'Stagger'
			for i in range(len(bp_idx)):
				if(len(self.base_pairs[bp_idx[i]].stagger)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_pairs[bp_idx[i]].stagger[midx])
					append = True

		if(parameter=='Buckle'):
			key = 'Buckle'
			for i in range(len(bp_idx)):
				if(len(self.base_pairs[bp_idx[i]].buckle)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_pairs[bp_idx[i]].buckle[midx])
					append = True

		if(parameter=='Propeller'):
			key = 'Propeller'
			for i in range(len(bp_idx)):
				if(len(self.base_pairs[bp_idx[i]].propeller)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_pairs[bp_idx[i]].propeller[midx])
					append = True

		if(parameter=='Opening'):
			key = 'Opening'
			for i in range(len(bp_idx)):
				if(len(self.base_pairs[bp_idx[i]].opening)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_pairs[bp_idx[i]].opening[midx])
					append = True
		
		if(parameter=='Radius S-1'):
			key = 'Radius S-1'
			for i in range(len(bp_idx)):
				if(len(self.base_pairs[bp_idx[i]].radS1)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_pairs[bp_idx[i]].radS1[midx])
					append = True


		if(parameter=='Radius S-2'):
			key = 'Radius S-2'
			for i in range(len(bp_idx)):
				if(len(self.base_pairs[bp_idx[i]].radS2)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_pairs[bp_idx[i]].radS2[midx])
					append = True

		#Base step parameters
		if(parameter=='Shift'):
			key = 'Shift'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].shift)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].shift[midx])
					append = True

		if(parameter=='Slide'):
			key = 'Slide'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].slide)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].slide[midx])
					append = True

		if(parameter=='Rise'):
			key = 'Rise'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].rise)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].rise[midx])
					append = True

		if(parameter=='Tilt'):
			key = 'Tilt'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].tilt)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].tilt[midx])
					append = True

		if(parameter=='Roll'):
			key = 'Roll'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].roll)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].roll[midx])
					append = True

		if(parameter=='Twist'):
			key = 'Twist'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].twist)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].twist[midx])
					append = True

		#Base step helical parameters
		if(parameter=='X-disp'):
			key = 'X-disp'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].x_disp)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].x_disp[midx])
					append = True

		if(parameter=='Y-disp'):
			key = 'Y-disp'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].y_disp)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].y_disp[midx])
					append = True

		if(parameter=='h-Rise'):
			key = 'h-Rise'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].h_rise)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].h_rise[midx])
					append = True

		if(parameter=='Inclination'):
			key = 'Inclination'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].inclination)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].inclination[midx])
					append = True

		if(parameter=='Tip'):
			key = 'Tip'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].tip)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].tip[midx])
					append = True

		if(parameter=='h-Twist'):
			key = 'h-Twist'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].h_twist)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].h_twist[midx])
					append = True
		
		# Helical axis related stuffs
		if(parameter=='Helical X-axis'):
			key = 'Helical X-axis'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].hel_Xaxis)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].hel_Xaxis[midx])
					append = True

		if(parameter=='Helical Y-axis'):
			key = 'Helical Y-axis'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].hel_Yaxis)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].hel_Yaxis[midx])
					append = True

		if(parameter=='Helical Z-axis'):
			key = 'Helical Z-axis'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].hel_Zaxis)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].hel_Zaxis[midx])
					append = True

		if(parameter=='Helical X-axis smooth'):
			key = 'Helical X-axis smooth'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].hel_Xaxis_smth)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].hel_Xaxis_smth[midx])
					append = True

		if(parameter=='Helical Y-axis smooth'):
			key = 'Helical Y-axis smooth'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].hel_Yaxis_smth)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].hel_Yaxis_smth[midx])
					append = True

		if(parameter=='Helical Z-axis smooth'):
			key = 'Helical Z-axis smooth'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].hel_Zaxis_smth)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].hel_Zaxis_smth[midx])
					append = True

		if(parameter=='Helical axis curvature'):
			key = 'Helical axis curvature'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].curvature)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].curvature[midx])
					append = True

		if(parameter=='Helical axis tangent'):
			key = 'Helical axis tangent'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].tangent)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].tangent[midx])
					append = True
		
		# Major and minor grooves
		if(parameter=='Major Groove'):
			key = 'Major Groove'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].major_pp)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].major_pp[midx])
					append = True
		
		if(parameter=='Major Groove Refined'):
			key = 'Major Groove Refined'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].major_refine)==0):
					empty = True
					key = 'Major Groove Refined'
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].major_refine[midx])
					append = True
		
		if(parameter=='Minor Groove'):
			key = 'Minor Groove'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].minor_pp)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].minor_pp[midx])
					append = True
		
		if(parameter=='Minor Groove Refined'):
			key = 'Minor Groove Refined'
			for i in range(len(bp_idx)):
				if(len(self.base_steps[bp_idx[i]].minor_refine)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_steps[bp_idx[i]].minor_refine[midx])
					append = True
		
		# Backbone dihedrals
		if(parameter=='alpha S-1'):
			key = 'alpha S-1'
			for i in range(len(bp_idx)):
				if(len(self.base_pairs[bp_idx[i]].alpha_s1)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_pairs[bp_idx[i]].alpha_s1[midx])
					append = True
		
		if(parameter=='beta S-1'):
			key = 'beta S-1'
			for i in range(len(bp_idx)):
				if(len(self.base_pairs[bp_idx[i]].beta_s1)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_pairs[bp_idx[i]].beta_s1[midx])
					append = True
		
		if(parameter=='gamma S-1'):
			key = 'gamma S-1'
			for i in range(len(bp_idx)):
				if(len(self.base_pairs[bp_idx[i]].gamma_s1)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_pairs[bp_idx[i]].gamma_s1[midx])
					append = True
		
		if(parameter=='delta S-1'):
			key = 'delta S-1'
			for i in range(len(bp_idx)):
				if(len(self.base_pairs[bp_idx[i]].delta_s1)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_pairs[bp_idx[i]].delta_s1[midx])
					append = True
		
		if(parameter=='epsilon S-1'):
			key = 'epsilon S-1'
			for i in range(len(bp_idx)):
				if(len(self.base_pairs[bp_idx[i]].epsilon_s1)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_pairs[bp_idx[i]].epsilon_s1[midx])
					append = True
		
		if(parameter=='zeta S-1'):
			key = 'zeta S-1'
			for i in range(len(bp_idx)):
				if(len(self.base_pairs[bp_idx[i]].zeta_s1)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_pairs[bp_idx[i]].zeta_s1[midx])
					append = True
		
		if(parameter=='chi S-1'):
			key = 'chi S-1'
			for i in range(len(bp_idx)):
				if(len(self.base_pairs[bp_idx[i]].chi_s1)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_pairs[bp_idx[i]].chi_s1[midx])
					append = True
		
		if(parameter=='alpha S-2'):
			key = 'alpha S-2'
			for i in range(len(bp_idx)):
				if(len(self.base_pairs[bp_idx[i]].alpha_s2)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_pairs[bp_idx[i]].alpha_s2[midx])
					append = True
		
		if(parameter=='beta S-2'):
			key = 'beta S-2'
			for i in range(len(bp_idx)):
				if(len(self.base_pairs[bp_idx[i]].beta_s2)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_pairs[bp_idx[i]].beta_s2[midx])
					append = True
		
		if(parameter=='gamma S-2'):
			key = 'gamma S-2'
			for i in range(len(bp_idx)):
				if(len(self.base_pairs[bp_idx[i]].gamma_s2)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_pairs[bp_idx[i]].gamma_s2[midx])
					append = True
		
		if(parameter=='delta S-2'):
			key = 'delta S-2'
			for i in range(len(bp_idx)):
				if(len(self.base_pairs[bp_idx[i]].delta_s2)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_pairs[bp_idx[i]].delta_s2[midx])
					append = True
		
		if(parameter=='epsilon S-2'):
			key = 'epsilon S-2'
			for i in range(len(bp_idx)):
				if(len(self.base_pairs[bp_idx[i]].epsilon_s2)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_pairs[bp_idx[i]].epsilon_s2[midx])
					append = True
		
		if(parameter=='zeta S-2'):
			key = 'zeta S-2'
			for i in range(len(bp_idx)):
				if(len(self.base_pairs[bp_idx[i]].zeta_s2)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_pairs[bp_idx[i]].zeta_s2[midx])
					append = True
		
		if(parameter=='chi S-2'):
			key = 'chi S-2'
			for i in range(len(bp_idx)):
				if(len(self.base_pairs[bp_idx[i]].chi_s2)==0):
					empty = True
					idx = bp_idx[i]
					break
				else:
					data.append(self.base_pairs[bp_idx[i]].chi_s2[midx])
					append = True
		
		if(empty):
			raise ValueError('ERROR: The parameter \"{0}\" for base pair/step \"{1}\" is not set/loaded.\n' .format(key,idx+1))
			exit(1)

		if(key=='dummy'):
			raise ValueError('ERROR: Incorrect parameter keyword: \"{0}\" .\n' .format(parameter))
			exit(1)
		
		return data, bp_idx

	def time_vs_parameter(self, parameter, bp, merge=False, merge_method='mean', masked=False):
		"""To get the parameter of either a specfic base-pair/step or a DNA segment as a function of time.
		
		**Arguments:**

			* ``parameter (string)``: Name of a base-pair or base-step or helical parameter
									For details about accepted keywords, see ``parameter`` in the method :meth:`dnaMD.DNA.get_parameters`.
			* ``bp (1D list) or (1D array)``: base-pairs to analyze
				Example: ::

							bp = [6]                                  # merge = False
							bp = [4,15]                               # merge = True
					
			* ``merge (bool)``: ``Dfault=False``: As shown above, if ``True``, bp should a list of range otherwise a list of single value. If ``bp = True``, the parameter for the respective DNA segment could be merged or calculated by ``merge_method``.

			* ``merge_method  (string)``: Method to calculate the parameter of a DNA segment from local parameters of all base-pairs/steps that are between the range given through ``bp``. 
				Currently accepted keywords are as follows:

					* ``merge_method = mean``: Average of local parameters
					* ``merge_method = sum``: Sum of local parameters


			* ``masked (bool)``: ``Dfault=False``: To skip specific frames/snapshots. dnaMD.DNA.mask array should be set to use this functionality. This array contains boolean (either ``True`` or ``False``) value for each frame to mask the frames. Presently, mask array is automatically generated during :meth:`dnaMD.DNA.generate_smooth_axis` method to skip those frames where 3D fitting curve was not successfull within the given critera.
           

        **Returns:**
            * ``time  (1D array)``: array containing time of length number of frames
            * ``value (1D array)``: array containing parameter values of length number of frames
		
		"""
		if not (isinstance(bp,list) or isinstance(bp,np.ndarray)):
			raise AssertionError("type %s is not list or np.ndarray" % type(bp))

		if (len(bp)>1) and (merge==False):
			raise AssertionError("bp %s contains more than two values, whereas merge=False. Use either one value in bp or merge=True" % bp)
			exit(1)

		if len(bp)==1:
			merge = False

		if (merge==True) and not ((merge_method == 'mean') or (merge_method == 'sum')):
			raise AssertionError("merge method %s is not available." % merge_method)
			exit(1)
		
		# Masking values according to mask array
		midx = []
		if masked and len(self.mask)==0:
			raise ValueError("mask array is not set. mask array is set within generate_smooth_axis() \n")
		for i in range(len(self.time)):
			if masked:
				if	self.mask[i] == False:
					midx.append(i)
			else:
				midx.append(i)


		if len(bp)==1:
			param_value, bp_idx = self.get_parameters(parameter,bp, bp_range=False, masked=masked)
		else:
			param_value, bp_idx = self.get_parameters(parameter,bp, bp_range=True, masked=masked)

		if (merge==True) and (merge_method=='mean'):
			return self.time, np.mean(param_value, axis=0)

		elif (merge==True) and (merge_method=='sum'):
			return self.time[midx], np.sum(param_value, axis=0)

		else:
			return self.time[midx], param_value[0]



	def parameter_distribution(self, parameter, bp, bins=30, merge=False, merge_method='mean', masked=False):
		"""To get the parameter distribution of either a specfic base-pair/step or a DNA segment over the MD simulation.
		
		**Arguments:**

			* ``parameter (string)``: Name of a base-pair or base-step or helical parameter
									For details about accepted keywords, see ``parameter`` in the method :meth:`dnaMD.DNA.get_parameters`.
					
			* ``bp (1D list) or (1D array)``: base-pairs to analyze
				Example: ::

							bp = [6]                                  # merge = False
							bp = [4,15]                               # merge = True
		
			* ``bins  (int)``: Number of bins to calculate histogram
			
			* ``merge (bool)``: ``Dfault=False``: As shown above, if ``True``, bp should a list of range otherwise a list of single value. If ``bp = True``, the parameter for the respective DNA segment could be merged or calculated by ``merge_method``.

			* ``merge_method  (string)``: Method to calculate the parameter of a DNA segment from local parameters of all base-pairs/steps that are between the range given through ``bp``. 
				Currently accepted keywords are as follows:

					* ``merge_method = mean``: Average of local parameters
					* ``merge_method = sum``: Sum of local parameters
           

			* ``masked (bool)``: ``Dfault=False``: To skip specific frames/snapshots. dnaMD.DNA.mask array should be set to use this functionality. This array contains boolean (either ``True`` or ``False``) value for each frame to mask the frames. Presently, mask array is automatically generated during :meth:`dnaMD.DNA.generate_smooth_axis` method to skip those frames where 3D fitting curve was not successfull within the given critera.


        **Returns:**
            * ``values   (1D array)``: array containing parameter values
            * ``density  (1D array)``: array containing density for respective parameter values
		
		"""
		if not (isinstance(bp,list) or isinstance(bp,np.ndarray)):
			raise AssertionError("type %s is not list or np.ndarray" % type(bp))

		if (len(bp)>1) and (merge==False):
			raise AssertionError("bp %s contains more than two values, whereas merge=False. Use either one value in bp or merge=True" % bp)
			exit(1)

		if len(bp)==1:
			merge = False

		if (merge==True) and not ((merge_method == 'mean') or (merge_method == 'sum')):
			raise AssertionError("merge method %s is not available." % merge_method)
			exit(1)

		if len(bp)==1:
			param_value, bp_idx = self.get_parameters(parameter,bp, bp_range=False, masked=masked)
		else:
			param_value, bp_idx = self.get_parameters(parameter,bp, bp_range=True, masked=masked)

		if (merge==True) and (merge_method=='mean'):
			param_value = np.mean(param_value, axis=0)

		elif (merge==True) and (merge_method=='sum'):
			param_value = np.sum(param_value, axis=0)

		else:
			param_value = param_value[0]
		
		density, bin_edges = np.histogram(param_value, bins=bins, density=True)
		bin_width = bin_edges[1]-bin_edges[0]

		density = np.insert(density, 0, 0.0)
		density = np.append(density, 0.0)

		values = []
		for i in range(len(bin_edges)-1):
			values.append((bin_edges[i]+bin_edges[i+1])/2)
		
		bin_width
		values = np.asarray(values)
		values = np.append(values, values[-1]+bin_width)
		values = np.insert(values, 0, values[0]-bin_width)

		return np.array(values), density

	def set_base_pair_parameters(self, filename, bp, parameters=[1,2,3,4,5,6], bp_range=True):
		"""	To read and store basepairs parameters (shear, stretch, stagger, buckle, propeller and opening) from an input file.
		
		**Arguments:**
		
			* ``filename (string)``: Input file, which is generated from do_x3dna. e.g. L-BP_g.dat
		    
			* ``bp (1D list) or (1D array)``: base-pairs to analyze
				Example: ::

							bp = [6]                                # bp_range = False
							bp = [4,15]                             # bp_range = True
							bp = range(4,15)                        # bp_range = False
							bp = np.arange(4,15)                    # bp_range = False
							bp = [2,5,6,7,9,12,18]                  # bp_range = False
							
			
			* ``parameters (1D list)``: List of numbers corrosponding to base-pairs parameters as follows:
						
						* ``shear      ->  1``
						* ``stretch    ->  2``
						* ``stagger    ->  3``
						* ``buckle     ->  4``
						* ``propeller  ->  5``
						* ``opening    ->  6``
				
				Example:
				
					*For shear, buckle, and propeller:*
					
								``parameters = [1,4,6]``
		            
					*For stretch, stagger, buckle, and propeller:*
		                
								``parameters = range(2,6)``
		                
								``parameters = [2,3,4,5]``
		
			* ``bp_range (bool)``: ``Dfault=True``: As shown above, if ``True``, bp is taken as a range otherwise list or numpy array

		**Returns:**
				``None``

		"""
		if not (isinstance(bp,list) or isinstance(bp,np.ndarray)):
			raise AssertionError("type %s is not list or np.ndarray" % type(bp))
		if not (isinstance(parameters,list) or isinstance(parameters, np.ndarray)):
			raise AssertionError("type %s is not list or np.ndarray" % type(parameters))

		for parameter in parameters:
			if(parameter>6):
				print '\nWARNING: number of requested parameters exceeded to {0} as contrast to six !!\n' .format(parameter)
				print 'Setting number of parameters from one to six\n\n'
				parameters = [1, 2, 3, 4, 5, 6]
				break

		data, time = read_param_file(filename,parameters,bp,bp_range)
		
		if(len(self.time)==0):
			self.time = np.array(time)
		else:
			if(len(time)!=len(self.time)):
				raise AssertionError("\nTime or number of frame mismatch in input files.\n Exiting...\n")
		
		bp_idx, param_idx = get_idx_of_bp_parameters(bp,parameters,bp_range)	
		
		for i in range(len(data)):
			for j in range(len(data[i])):
				if(0==param_idx[j]):
					self.base_pairs[bp_idx[i]].shear = data[i][j]
				if(1==param_idx[j]):
					self.base_pairs[bp_idx[i]].stretch = data[i][j]
				if(2==param_idx[j]):
					self.base_pairs[bp_idx[i]].stagger = data[i][j]
				if(3==param_idx[j]):
					self.base_pairs[bp_idx[i]].buckle = data[i][j]
				if(4==param_idx[j]):
					self.base_pairs[bp_idx[i]].propeller = data[i][j]
				if(5==param_idx[j]):
					self.base_pairs[bp_idx[i]].opening = data[i][j]


	def set_major_minor_groove(self, filename, bp_step, parameters=[1,2,3,4], step_range=True):
		"""	To read and store Major and Minor grooves from an input file.

			* Minor groove : direct P-P distance
			* Minor Grrove Refined : refined P-P distance which take into account the directions of the sugar-phosphate backbones
			* Major groove : direct P-P distance
			* Major Grrove Refined : refined P-P distance which take into account the directions of the sugar-phosphate backbones
			
			.. warning::

				* The major and minor grooves (direct P-P) cannot be calculated for first and last two base-steps
				* The major and minor grooves (refined P-P) cannot be calculated for first and last three base-steps

		**Arguments:**
		
			* ``filename (string)``: Input file, which is generated from do_x3dna. e.g. L-BP_g.dat
		    
			* ``bp_step (1D list) or (1D array)``: bases-steps to analyze
				Example: ::

							bp_step = [6]                                # step_range = False
							bp_step = [4,15]                             # step_range = True
							bp_step = range(4,15)                        # step_range = False
							bp_step = np.arange(4,15)                    # step_range = False
							bp_step = [2,5,6,7,9,12,18]                  # step_range = False
							
			
			* ``parameters (1D list)``: List of numbers corrosponding to base-pairs parameters as follows:
						
						* ``Minor Groove          ->  1``
						* ``Minor Grrove Refined  ->  2``
						* ``Major Groove          ->  3``
						* ``Major Grrove Refined  ->  4``
				
				Example:
				
					*For minor (refined) and major (refined) grooves:*
					
								``parameters = [2,4]``
		            
		
			* ``step_range (bool)``: ``Dfault=True``: As shown above, if ``True``, bp_step is taken as a range otherwise list or numpy array

		**Returns:**
				``None``

		"""
		if not (isinstance(bp_step,list) or isinstance(bp_step,np.ndarray)):
			raise AssertionError("type %s is not list or np.ndarray" % type(bp_step))
		if not (isinstance(parameters,list) or isinstance(parameters, np.ndarray)):
			raise AssertionError("type %s is not list or np.ndarray" % type(parameters))

		for parameter in parameters:
			if(parameter>4):
				print '\nWARNING: number of requested parameters exceeded to {0} as contrast to four !!\n' .format(parameter)
				print 'Setting number of parameters from one to six\n\n'
				parameters = [1, 2, 3, 4]
				break

		data, time = read_param_file(filename, parameters, bp_step, step_range, word=True)
		
		if(len(self.time)==0):
			self.time = np.array(time)
		else:
			if(len(time)!=len(self.time)):
				raise AssertionError("\nTime or number of frame mismatch in input files.\n Exiting...\n")

		
		bp_idx, param_idx = get_idx_of_bp_parameters(bp_step, parameters, step_range)	
		
		for i in range(len(data)):
			for j in range(len(data[i])):
				if (0==param_idx[j]) and (data[i][j][0] != None):
					self.base_steps[bp_idx[i]].minor_pp = data[i][j]
						
				if (1==param_idx[j]) and (data[i][j][0] != None):
					self.base_steps[bp_idx[i]].minor_refine = data[i][j]
				
				if (2==param_idx[j]) and (data[i][j][0] != None):
					self.base_steps[bp_idx[i]].major_pp = data[i][j]
				
				if (3==param_idx[j]) and (data[i][j][0] != None):
					self.base_steps[bp_idx[i]].major_refine = data[i][j]

	def set_backbone_dihedrals(self, filename, bp, parameters=range(1,15), bp_range=True):
		"""	To read and store backbone dihedrals (alpha, beta, gamma, delta, epsilon and zeta) and chi dihedral of both strands from an input file.

		.. note::

			* alpha:   O3'(i-1)-P-O5'-C5'
			* beta:    P-O5'-C5'-C4'
			* gamma:   O5'-C5'-C4'-C3'
			* delta:   C5'-C4'-C3'-O3'
			* epsilon: C4'-C3'-O3'-P(i+1)
			* zeta:    C3'-O3'-P(i+1)-O5'(i+1)
			* chi for pyrimidines(Y): O4'-C1'-N1-C2
			* chi for purines(R): O4'-C1'-N9-C4

		
		**Arguments:**
		
			* ``filename (string)``: Input file, which is generated from do_x3dna. e.g. L-BP_g.dat
		    
			* ``bp (1D list) or (1D array)``: base-pairs to analyze
				Example: ::

							bp = [6]                                # bp_range = False
							bp = [4,15]                             # bp_range = True
							bp = range(4,15)                        # bp_range = False
							bp = np.arange(4,15)                    # bp_range = False
							bp = [2,5,6,7,9,12,18]                  # bp_range = False
							
			
			* ``parameters (1D list)``: List of numbers corrosponding to base-pairs parameters as follows:
						
						* ``Alpha Strand I    ->   1``
						* ``Beta Strand I     ->   2``
						* ``Gamma Strand I    ->   3``
						* ``Delta Strand I    ->   4``
						* ``Epsilon Strand I  ->   5``
						* ``Zeta Strand I     ->   6``
						* ``Chi Strand I      ->   7``
						* ``Alpha Strand II   ->   8``
						* ``Beta Strand II    ->   9``
						* ``Gamma Strand II   ->  10``
						* ``Delta Strand II   ->  11``
						* ``Epsilon Strand II ->  12``
						* ``Zeta Strand II    ->  13``
						* ``Chi Strand II     ->  14``
				
				Example:
				
					*For alpha, delta, and zeta of strand I:*
					
								``parameters = [1,4,6]``
		            
					*For beta, gamma, delta, and epsilon of strand II:*
		                
								``parameters = range(9,13)``
		                
								``parameters = [9, 10, 11, 12]``
		
			* ``bp_range (bool)``: ``Dfault=True``: As shown above, if ``True``, bp is taken as a range otherwise list or numpy array

		**Returns:**
				``None``

		"""
		if not (isinstance(bp,list) or isinstance(bp,np.ndarray)):
			raise AssertionError("type %s is not list or np.ndarray" % type(bp))
		if not (isinstance(parameters,list) or isinstance(parameters, np.ndarray)):
			raise AssertionError("type %s is not list or np.ndarray" % type(parameters))

		for parameter in parameters:
			if(parameter>14):
				print '\nWARNING: number of requested parameters exceeded to {0} as contrast to six !!\n' .format(parameter)
				print 'Setting number of parameters from one to six\n\n'
				parameters = range(1,15)
				break

		data, time = read_param_file(filename,parameters,bp,bp_range,word=True)
		
		if(len(self.time)==0):
			self.time = np.array(time)
		else:
			if(len(time)!=len(self.time)):
				raise AssertionError("\nTime or number of frame mismatch in input files.\n Exiting...\n")
		
		bp_idx, param_idx = get_idx_of_bp_parameters(bp,parameters,bp_range)	
		
		for i in range(len(data)):
			for j in range(len(data[i])):
				if (0==param_idx[j]) and (data[i][j][0] != None):
					self.base_pairs[bp_idx[i]].alpha_s1 = data[i][j]
				if (1==param_idx[j]) and (data[i][j][0] != None):
					self.base_pairs[bp_idx[i]].beta_s1 = data[i][j]
				if (2==param_idx[j]) and (data[i][j][0] != None):
					self.base_pairs[bp_idx[i]].gamma_s1 = data[i][j]
				if (3==param_idx[j]) and (data[i][j][0] != None):
					self.base_pairs[bp_idx[i]].delta_s1 = data[i][j]
				if (4==param_idx[j]) and (data[i][j][0] != None):
					self.base_pairs[bp_idx[i]].epsilon_s1 = data[i][j]
				if (5==param_idx[j]) and (data[i][j][0] != None):
					self.base_pairs[bp_idx[i]].zeta_s1 = data[i][j]
				if (6==param_idx[j]) and (data[i][j][0] != None):
					self.base_pairs[bp_idx[i]].chi_s1 = data[i][j]
				if (7==param_idx[j]) and (data[i][j][0] != None):
					self.base_pairs[bp_idx[i]].alpha_s2 = data[i][j]
				if (8==param_idx[j]) and (data[i][j][0] != None):
					self.base_pairs[bp_idx[i]].beta_s2 = data[i][j]
				if (9==param_idx[j]) and (data[i][j][0] != None):
					self.base_pairs[bp_idx[i]].gamma_s2 = data[i][j]
				if (10==param_idx[j]) and (data[i][j][0] != None):
					self.base_pairs[bp_idx[i]].delta_s2 = data[i][j]
				if (11==param_idx[j]) and (data[i][j][0] != None):
					self.base_pairs[bp_idx[i]].epsilon_s2 = data[i][j]
				if (12==param_idx[j]) and (data[i][j][0] != None):
					self.base_pairs[bp_idx[i]].zeta_s2 = data[i][j]
				if (13==param_idx[j]) and (data[i][j][0] != None):
					self.base_pairs[bp_idx[i]].chi_s2 = data[i][j]



	def set_helical_radius(self, filename, bp, atomname='P', full=False, bp_range=True):
		"""	To read and set local helical radius of both strand
		
		**Arguments:**
			
			* ``filename (string)``: Input file, which is generated from do_x3dna. e.g. HelixRad_g.dat
			
			* ``bp (1D list) or (1D array)``: base-pairs to analyze
				Example: ::

							bp = [6]                                # bp_range = False
							bp = [4,15]                             # bp_range = True
							bp = range(4,15)                        # bp_range = False
							bp = np.arange(4,15)                    # bp_range = False
							bp = [2,5,6,7,9,12,18]                  # bp_range = False
			

			* ``atomname (string)``: list of atom names to consider for the DNA helix (accepted keywords: ``P, O4*, O4', C1* and C1``)
			
			* ``full (bool)``: To calculate full helical radius. Overrides atomname option and uses atom 'P', subsequently added 1 A to the radius calculated by 3DNA package
			
			* ``bp_range (bool)``: Shown above. if True, bp should be a range otherwise list or numpy array

		**Returns:** 
				``None``
		"""
		if not (isinstance(bp,list) or isinstance(bp,np.ndarray)):
			raise AssertionError("type %s is not list or np.ndarray" % type(bp))

		if not ( (atomname=='P') or (atomname=='O4*') or (atomname=='C1*') or (atomname=='O4\'') or (atomname=='C1\'')):
			print '\n This atomname {0} is not implemented... Exiting\n' .format(atomname)

		
		parameter = []
		if (atomname=='P') or (full):
			parameter = [1,4]
		
		if (atomname=='O4*' or atomname=='O4\'') and (not full):
			parameter = [2,5]

		if (atomname=='C1*' or atomname=='C1\'') and (not full):
			parameter = [3,6]

		data, time = read_param_file(filename,[1,2,3,4,5,6], bp, bp_range)
		
		if(len(self.time)==0):
			self.time = np.array(time)
		else:
			if(len(time)!=len(self.time)):
				raise AssertionError("\nTime or number of frame mismatch in input files.\n Exiting...\n")
		
		bp_idx, param_idx = get_idx_of_bp_parameters(bp,parameter, bp_range)
		
		if full:
			data = np.add(data, 1.0)
		
		for i in range(len(data)):
			if (atomname=='P') or (full):
				self.base_pairs[bp_idx[i]].radS1 = data[i][0]
				self.base_pairs[bp_idx[i]].radS2 = data[i][3]

			if (atomname=='O4*' or atomname=='O4\''):
				self.base_pairs[bp_idx[i]].radS1 = data[i][1]
				self.base_pairs[bp_idx[i]].radS2 = data[i][4]

			if (atomname=='C1*' or atomname=='C1\''):
				self.base_pairs[bp_idx[i]].radS1 = data[i][2]
				self.base_pairs[bp_idx[i]].radS2 = data[i][5]

	def set_base_step_parameters(self, filename, bp_step, parameters=[1,2,3,4,5,6], step_range=True,helical=False):
		"""	To read and store base-step (Shift, Slide, Rise, Tilt, Roll and Twist) and helical base-step (X-disp, Y-disp, h-Rise, Inclination, Tip and h-Twist) parameters from an input file
		
		**Arguments:**
		
			* ``filename (string)``: Input file, which is generated from do_x3dna. e.g. ``L-BPS_g.dat`` or ``L-BPH_g.dat``
		    
			* ``bp_step (1D list) or (1D array)``: base-steps to analyze
				Example: ::

							bp_step = [6]                                # step_range = False
							bp_step = [4,15]                             # step_range = True
							bp_step = range(4,15)                        # step_range = False
							bp_step = np.arange(4,15)                    # step_range = False
							bp_step = [2,5,6,7,9,12,18]                  # step_range = False
							
			
			* ``parameters (1D list)``: List of numbers corrosponding to base-steps parameters as follows:
		            
					If ``helical = False``:
							* ``Shift       ->  1``
							* ``Slide       ->  2``
							* ``Rise        ->  3``
							* ``Tilt        ->  4``
							* ``Roll        ->  5``
							* ``Twist       ->  6``

					If ``helical = True``:
							* ``X-disp      ->  1``
							* ``Y-disp      ->  2``
							* ``h-Rise      ->  3``
							* ``Inclination ->  4``
							* ``Tip         ->  5``
							* ``h-Twist     ->  6``
					
					Example:
						
						*For Shift, Tilt and Twist*:
						
								``parameters = [1,4,6]``
								
						*For Slide, Rise, Tilt and Roll*:
						
								``parameters = range(2,6)``
								
								``parameters = [2,3,4,5]``
			
			* ``step_range (bool)``: ``Dfault=True``: As shown above, if ``True``, bp_step is taken as a range otherwise list or numpy array

		**Returns:**
				``None``
		"""
		if not (isinstance(bp_step,list) or isinstance(bp_step,np.ndarray)):
			raise AssertionError("type %s is not list or np.ndarray" % type(bp_step))
		if not (isinstance(parameters,list) or isinstance(parameters, np.ndarray)):
			raise AssertionError("type %s is not list or np.ndarray" % type(parameters))

		for parameter in parameters:
			if(parameter>6):
				print '\nWARNING: number of requested parameters exceeded to {0} as contrast to six !!\n' .format(parameter)
				print 'Setting number of parameters from one to six\n\n'
				parameters=range(1,6)
				break

		data, time = read_param_file(filename,parameters,bp_step,step_range)
		
		if(len(self.time)==0):
			self.time = np.array(time)
		else:
			if(len(time)!=len(self.time)):
				raise AssertionError("\nTime or number of frame mismatch in input files.\n Exiting...\n")
		
		bp_idx, param_idx = get_idx_of_bp_parameters(bp_step,parameters,step_range)	
		
		for i in range(len(data)):
			for j in range(len(data[i])):
				if(not helical):
					if(0==param_idx[j]):
						self.base_steps[bp_idx[i]].shift = data[i][j]
					if(1==param_idx[j]):
						self.base_steps[bp_idx[i]].slide = data[i][j]
					if(2==param_idx[j]):
						self.base_steps[bp_idx[i]].rise = data[i][j]
					if(3==param_idx[j]):
						self.base_steps[bp_idx[i]].tilt = data[i][j]
					if(4==param_idx[j]):
						self.base_steps[bp_idx[i]].roll = data[i][j]
					if(5==param_idx[j]):
						self.base_steps[bp_idx[i]].twist = data[i][j]
				if(helical):
					if(0==param_idx[j]):
						self.base_steps[bp_idx[i]].x_disp = data[i][j]
					if(1==param_idx[j]):
						self.base_steps[bp_idx[i]].y_disp = data[i][j]
					if(2==param_idx[j]):
						self.base_steps[bp_idx[i]].h_rise = data[i][j]
					if(3==param_idx[j]):
						self.base_steps[bp_idx[i]].inclination = data[i][j]
					if(4==param_idx[j]):
						self.base_steps[bp_idx[i]].tip = data[i][j]
					if(5==param_idx[j]):
						self.base_steps[bp_idx[i]].h_twist = data[i][j]
		
	def get_mean_error(self,bp,parameter,err_type='std',bp_range=True, merge_bp=1, merge_method='mean', masked=False):
		"""To calculate average and error of the given parameter for the gieven set of base-pairs/steps
		
		**Arguments:**
			* ``bp (1D list) or (1D array)``: base-pairs or base-steps to analyze
				Example: ::
					
					bp = [6]                                # bp_range = False
					bp = [4,15]                             # bp_range = True
					bp = range(4,15)                        # bp_range = False
					bp = np.arange(4,15)                    # bp_range = False
			
			* ``parameter (string)``: Name of a base-pair or base-step or helical parameter
									For details about accepted keywords, see ``parameter`` in the method :meth:`DNA.get_parameters`.
			
			* ``error     (string)``:  Method of error estimation.
				Currently accepted method as follows:
					* ``error = 'std'``   : Standard Deviation
					* ``error = 'acf'``   : Standard error using autocorrelation time (requires: g_analyze)
					* ``error = 'block'`` : Standard error using block averaging method (requires: g_analyze)

					.. warning::
						to calculate errors by using ``error = 'acf'`` or ``error = 'block'`` , GROMACS tool ``g_analyze`` should be present in ``$PATH``.
			
			* ``bp_range (bool)``: Shown above. if True, bp should be a range otherwise list or numpy array
			
			* ``merge_bp (integer)``: Number of base-pairs or steps to merge for creating the small DNA segments
			
			* ``merge_method  (string)``: Method to calculate the parameter of a DNA segment from local parameters of all base-pairs/steps that are between the range given through ``bp``. 
				Currently accepted keywords are as follows:

					* ``merge_method = mean``: Average of local parameters
					* ``merge_method = sum``: Sum of local parameters
            
			* ``masked (bool)``: ``Dfault=False``: To skip specific frames/snapshots. dnaMD.DNA.mask array should be set to use this functionality. This array contains boolean (either ``True`` or ``False``) value for each frame to mask the frames. Presently, mask array is automatically generated during :meth:`dnaMD.DNA.generate_smooth_axis` method to skip those frames where 3D fitting curve was not successfull within the given critera.
		

		**Returns:**
				* ``basepairs/steps (1D array)``: Number of base pair-steps. If ``merge_bp>1``, middle number will be returned.
				* ``avg. parameter values (1D array)``: average values of the parameter
				* ``error (1D array)``: Error values for corresponding average values
		
		"""

		merge=False

		if(not bp_range) and (merge_bp>1):
			print ("\nERROR: Merging of base pairs/steps only possible with given base pair/steps range\n")
			exit(1)

		if(bp_range) and (merge_bp>1):
			merge = True
		
		if (merge==True) and not ((merge_method == 'mean') or (merge_method == 'sum')):
			raise AssertionError("merge method %s is not available." % merge_method)
			exit(1)
		
		data, bp_idx = self.get_parameters(parameter,bp,bp_range, masked=masked)

		bp_number = np.add(bp_idx,1)
		data = np.array(data)

		# Merging base pairs/step data for the given parameters
		merge_data = []
		mid_bin = []
		if(merge):
			i=0
			while(1):
				start = i
				if((i + merge_bp)>=len(bp_idx)):
					end = i+(len(bp_idx)-i)
				else:
					end   = i + merge_bp	
				
				if(start < end):
					mid_bin.append(int(start+(end-start)/2) + bp_number[0])
					if (merge_method == 'mean'):
						merge_data.append(np.mean(data[start:end],axis=0))
					if (merge_method == 'sum'):
						merge_data.append(np.sum(data[start:end],axis=0))
				
				if(i >= len(bp_idx)):
					break
				
				i += merge_bp

			merge_data = np.array(merge_data)
	
		
		if (err_type == 'std'):
			if(merge):
				error = np.std(merge_data,axis=1)
			else:
				error = np.std(data,axis=1)

		if (err_type =='acf' ) or (err_type =='block' ):
			if(merge):
				error = get_error(self.time, merge_data, len(mid_bin), err_type=err_type)
			else:
				error = get_error(self.time, data, len(bp_idx), err_type=err_type)

		if(merge):
			return mid_bin, np.mean(merge_data,axis=1), error
		else:
			return bp_number, np.mean(data,axis=1), error

	def set_helical_axis(self, filename, step_range=False,step=None):
		"""
		To read and set local helical-axis postions from an input file.
	
		**Arguments:**
		
			* ``filename (string)``: Input file, which is generated from do_x3dna. e.g. HelAxis_g.dat
			
			* ``step_range (bool)``:
					* ``step_range = True`` : read axis coordinates of base-steps for the given range of base-steps
					* ``step_range = False``: read axis coordinates of all base-steps
					
			* ``step (list)``: list containing lower and higher limit of base-steps range
							* This option only works with ``step_range=True``.
							* This list should not contain more than two number. 
							* First number should be less than second number.
							
					Example:
						
						*For base-step 4 to 15*:
							``step = [4,15]         # step_range = True``
							
		**Returns:**
					``None``

		"""

		if (step_range):
			if not isinstance(step,list):
				raise AssertionError("type %s is not list" % type(step))
			if (len(step)>2):
				print ("ERROR: Range for helical axis should be list of two numbers, e.g. step=[1, 20] \n")
				exit(1)
		
		if (step_range) and (step == None):
			raise ValueError("See, documentation for step  and step_range usage!!!")
			

		if (step_range):
			if (len(step) != 2):
				raise ValueError("See, documentation for step usage!!!")
			
			if step[0] > step[1]:
				raise ValueError("See, documentation for step usage!!!")
			data, time = read_param_file(filename, [1,2,3], step, True)
		else:
			data, time = read_param_file(filename, [1,2,3], [1,self.num_step], True)
		
		if(len(self.time)==0):
			self.time = np.array(time)
		else:
			if(len(time)!=len(self.time)):
				raise AssertionError("\nTime or number of frame mismatch in input files.\n Exiting...\n")

		if (step_range):
			bp_idx, param_idx = get_idx_of_bp_parameters(step,[],True)
		else:
			bp_idx, param_idx = get_idx_of_bp_parameters([1,self.num_step],[],True)

		for i in range(len(data)):
			for j in range(len(data[i])):
				if(0==j):
					self.base_steps[bp_idx[i]].hel_Xaxis = data[i][j]
				if(1==j):
					self.base_steps[bp_idx[i]].hel_Yaxis = data[i][j]
				if(2==j):
					self.base_steps[bp_idx[i]].hel_Zaxis = data[i][j]
	
	def generate_smooth_axis(self, step_range=False, step=None, smooth=500.0, spline=3, fill_point=6, cut_off_angle=20):
		"""	To smoothen the helical axis using spline interpolation.

		.. note::
			A 3D curve is fitted on local helical axis that are calculated using ``do_x3dna`` tool. Sometimes in few frames, fitting **may not** be accurate and produces artifact. To record these frames, dnaMD.DNA.mask array containing boolean values are generated. If value is ``True``, fitting might not be correct and vice versa. This array could be used in later analysis to skip/mask the frames containing inaccurate axis.

		.. warning::
			This function requires `SciPy package <http://www.scipy.org/>`_.
		
		**Arguments:**
			
			* ``step_range (bool)``:
					* ``step_range = True`` : Smoothen axis for the given range of base-steps
					* ``step_range = False``: Smoothen axis for entire DNA. If original helical-axis of any base-step will be found to be not available, error will be raised.
					
			* ``step (list)``: list containing lower and higher limit of base-steps range
							* This option only works with ``step_range=True``.
							* This list should not contain more than two number. 
							* First number should be less than second number.
							
					Example:
						
						*For base-step 4 to 15*:
							``step = [4,15]         # step_range = True``

			* ``smooth (float)``: A smoothing condition. For more details, see about ``s = None``, which is paased into  `scipy.interpolate.splprep() <http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.splprep.html#scipy.interpolate.splprep>`_ method. 
				.. warning ::

					* Lower value may lead to an artifact of local sharp kink in the smoothed axis. 
					* Higher value may lead to the calculation of wrong helical axis.


			* ``spline (int)``: Degree of spline. For more details, see about ``k = 3``, which is paased into  `scipy.interpolate.splprep() <http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.splprep.html#scipy.interpolate.splprep>`_ method.

			
			* ``fill_point (int)``: Number of intrapolated points between two adjacent helical-axis coordinates.
			
			
			* ``cut_off_angle (float)``: Cut-off bending angle to define sharp kink in fitted curve. If angle in fitted curve is larger than this cut-off, refitting will be performed after deleting few of the original helical axis positions. If after this deletions, bending angle will not reduce below cut-off angle, value of ``smooth`` will be increased by 100 and entire cycle of fitting-refitting will be performed. When, value of ``smooth`` increases to more than 10000 during this fitting-refitting cycles, fitting process will be stopped with a warning message.
			
			
		**Returns:**
					``None``

		"""
		
		if not scipy_imported:
			raise ImportError("SciPy package is not available. Please visit http://www.scipy.org/install.html for download and installation instructions.\n") 

		if (step_range) and (step == None):
			raise ValueError("See, documentation for step  and step_range usage!!!")
		
		bp_idx = []	
		
		if step_range:
			if (len(step) != 2):
				raise ValueError("See, documentation for step usage!!!")
			
			if step[0] > step[1]:
				raise ValueError("See, documentation for step usage!!!")
			RawX, bp_idx = self.get_parameters('Helical X-axis', step, bp_range=True)
			RawY, dummy = self.get_parameters('Helical Y-axis', step, bp_range=True)
			RawZ, dummy = self.get_parameters('Helical Z-axis', step, bp_range=True)
		else:
			RawX, bp_idx = self.get_parameters('Helical X-axis', [1, self.num_step], bp_range=True)
			RawY, dummy = self.get_parameters('Helical Y-axis', [1, self.num_step], bp_range=True)
			RawZ, dummy = self.get_parameters('Helical Z-axis', [1, self.num_step], bp_range=True)

		RawX = np.array(RawX).T
		RawY = np.array(RawY).T
		RawZ = np.array(RawZ).T

		smoothX, smoothY, smoothZ = [], [], []

		nframes = len(self.time)
		
		if len(self.mask)==0:
			self.mask = np.zeros(nframes, dtype=bool)

		for i in range(nframes):
			
			frame_number = i+1	
			if((frame_number<100) and (frame_number%10==0)) or ((frame_number<1000) and (frame_number%100==0)) or (frame_number%1000==0):
				sys.stdout.write("\rFitting spline curve on helcial axis of frame %d out of %d frames" % (frame_number, nframes))
				sys.stdout.flush()

			xsmooth, ysmooth, zsmooth, mask = fit_axis(bp_idx, frame_number, RawX[i], RawY[i], RawZ[i], smooth, spline, fill_point, cut_off_angle)
			self.mask[i] = mask

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
			self.base_steps[bp_idx[i]].hel_Xaxis_smth = smoothX[i]
			self.base_steps[bp_idx[i]].hel_Yaxis_smth = smoothY[i]
			self.base_steps[bp_idx[i]].hel_Zaxis_smth = smoothZ[i]

	def write_haxis_pdb(self, filename='helical_axis.pdb', step_range=False, step=None, write_smooth_axis=True, write_orig_axis=False, write_curv=False, scale_curv=1):
		"""
		To write trajectory of helcial-axis as a PDB format file. Helical axis could be original or smoothed. For smoothed axis, curvature could be written in B-factor field of PDB file.

		**Arguments:**
			
			* ``filename (string)``: Name of the output PDB format file.
			
			* ``step_range (bool)``:
					* ``step_range = True`` : Smoothen axis for the given range of base-steps
					* ``step_range = False``: Smoothen axis for entire DNA. If original helical-axis of any base-step will be found to be not available, error will be raised.
					
			* ``step (list)``: list containing lower and higher limit of base-steps range
							* This option only works with ``step_range=True``.
							* This list should not contain more than two number. 
							* First number should be less than second number.
							
					Example:
						
						*For base-step 4 to 15*:
							``step = [4,15]         # step_range = True``

			* ``write_smooth_axis (bool)``: Write coordinates of smoothed helical axis as chain A.
			
			* ``write_orig_axis (bool)``: Write coordinates of original helical axis (output from do_x3dna) as chain B.
			
			* ``write_curv (bool)``: Write curvature of smoothed helical axis in B-factor coloumn of PDB file.
			
			* ``scale_curv (int)``: Scaling of curvature. ``curvature * scale_curv`` is written in  B-factor coloumn of PDB file.
			
		**Returns:**

				None

		"""
		
		if (step_range) and (step == None):
			raise ValueError("See, documentation for step  and step_range usage!!!")
		
		if not write_orig_axis and not write_smooth_axis:
			raise ValueError("Nothing to write as both \"write_orig_axis=Flase\" and \"write_smooth_axis=False\" !!!")


		if step_range:
			if (len(step) != 2):
				raise ValueError("See, documentation for step usage!!!")
			
			if step[0] > step[1]:
				raise ValueError("See, documentation for step usage!!!")
			
			# Orignal helical axis
			if (write_orig_axis):
				RawX, bp_idx = self.get_parameters('Helical X-axis', step, bp_range=True)
				RawY, dummy = self.get_parameters('Helical Y-axis', step, bp_range=True)
				RawZ, dummy = self.get_parameters('Helical Z-axis', step, bp_range=True)
			
			# Smoothed helical axis
			if (write_smooth_axis):
				SmoothX, bp_idx = self.get_parameters('Helical X-axis smooth', step, bp_range=True)
				SmoothY, bp_idx = self.get_parameters('Helical Y-axis smooth', step, bp_range=True)
				SmoothZ, bp_idx = self.get_parameters('Helical Z-axis smooth', step, bp_range=True)

			# Helical axis curvature
			if (write_curv):
				curvature, bp_idx = self.get_parameters('Helical axis curvature', step, bp_range=True)
		
		else:

			# Orignal helical axis
			if (write_orig_axis):
				RawX, bp_idx = self.get_parameters('Helical X-axis', [1, self.num_step], bp_range=True)
				RawY, dummy = self.get_parameters('Helical Y-axis', [1, self.num_step], bp_range=True)
				RawZ, dummy = self.get_parameters('Helical Z-axis', [1, self.num_step], bp_range=True)
			
			# Smoothed helical axis
			if (write_smooth_axis):
				SmoothX, bp_idx = self.get_parameters('Helical X-axis smooth', [1, self.num_step], bp_range=True)
				SmoothY, bp_idx = self.get_parameters('Helical Y-axis smooth', [1, self.num_step], bp_range=True)
				SmoothZ, bp_idx = self.get_parameters('Helical Z-axis smooth', [1, self.num_step], bp_range=True)
			
			# Helical axis curvature
			if (write_curv):
				curvature, bp_idx = self.get_parameters('Helical axis curvature', [1, self.num_step], bp_range=True)
	
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
			f.write('%-6s    %4d\n' % ("MODEL",i+1))
			
			bfactor = 0.00

			if (write_smooth_axis):
				for j in range(len(SmoothX[i])):
					
					if (write_curv):
						bfactor = curvature[i][j]*scale_curv
					
					f.write('%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n' % ("ATOM",j+1,"CA"," ","AXS","A",j+1, " ", SmoothX[i][j], SmoothY[i][j], SmoothZ[i][j], 1.00, bfactor))

				for j in range(len(SmoothX[i])-1):
					f.write('CONECT %4d %4d\n' % (j+1, j+2))
				
				f.write("TER\n")
				
			if (write_orig_axis):
				atomstart = 0
				if (write_smooth_axis):
					atomstart = len(SmoothX[i])

				for j in range(len(RawX[i])):
					f.write('%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n' % ("ATOM",j+1+atomstart,"O"," ","AXS","B",j+1, " ", RawX[i][j],RawY[i][j],RawZ[i][j], 1.00, 0.00))

				for j in range(len(RawX[i])-1):
					f.write('CONECT %4d %4d\n' % (j+1+atomstart, j+2+atomstart))
				
				f.write("TER\n")

			f.write("ENDMDL\n")
	
		f.close()


	def calculate_curvature_tangent(self, step_range=False, step=None, store_tangent=False):
		"""	
		To calculate curvatures and tangent vectors along the helical axis. The curvature and tangent vectors are calculated using Frenet-Serret formula. The calculated values are stored in DNA object.
		
		**Arguments:**
			
			* ``step_range (bool)``:
					* ``step_range = True`` : Calculate curvature and tangent vectors for the given range of base-steps
					* ``step_range = False``: Calculate curvature and tangent vectors for entire DNA. If smoothed helical-axis of any base-step will be found to be not available, error will be raised.
					
			* ``step (list)``: list containing lower and higher limit of base-steps range
							* This option only works with ``step_range=True``.
							* This list should not contain more than two number. 
							* First number should be less than second number.
							
					Example:
						
						*For base-step 4 to 15*:
							``step = [4,15]         # step_range = True``

			
			* ``store_tangent (bool)``:
					* ``store_tangent = True`` : The calculated tangent vectors will be stored for later use.
					* ``store_tangent = False``:  The calculated tangent vectors will be discarded.
			

		**Returns:**
					``None``

		"""
		
		
		if not self.smooth_axis:
			raise ValueError("The helical axis is not smooth. At first, smooth the axis using generate_smooth_axis() method as described in http://rjdkmr.github.io/do_x3dna/apidoc.html#dnaMD.DNA.generate_smooth_axis.")

		if (step_range) and (step == None):
			raise ValueError("See, documentation for step  and step_range usage!!!")
		
		if step_range:
			if (len(step) != 2):
				raise ValueError("See, documentation for step usage!!!")
			
			if step[0] > step[1]:
				raise ValueError("See, documentation for step usage!!!")
			
			X, bp_idx = self.get_parameters('Helical X-axis smooth', step, bp_range=True)
			Y, bp_idx = self.get_parameters('Helical Y-axis smooth', step, bp_range=True)
			Z, bp_idx = self.get_parameters('Helical Z-axis smooth', step, bp_range=True)
		else:
			X, bp_idx = self.get_parameters('Helical X-axis smooth', [1, self.num_step], bp_range=True)
			Y, bp_idx = self.get_parameters('Helical Y-axis smooth', [1, self.num_step], bp_range=True)
			Z, bp_idx = self.get_parameters('Helical Z-axis smooth', [1, self.num_step], bp_range=True)

		X = np.asarray(X).T
		Y = np.asarray(Y).T
		Z = np.asarray(Z).T

		curvature, tangent = [], []
		
		for i in range(len(self.time)):
		
			#Curvature calculation
			xyz=np.vstack((X[i], Y[i], Z[i])).T
			T, N, B, k_temp, t_temp = frenet_serret(xyz)
			
			curvature.append(k_temp.flatten())
			
			if(store_tangent):
				tangent.append(T)
			

		curvature = np.asarray(curvature).T
		for i in range(len(bp_idx)):
			self.base_steps[bp_idx[i]].curvature = curvature[i]

		if(store_tangent):
			tangent = np.asarray(tangent)
			final_tan = []
			for i in range(len(tangent[0])):
				temp = []
				for j in range(len(tangent)):
					temp.append(tangent[j][i])
				final_tan.append(np.asarray(temp))
			
			for i in range(len(bp_idx)):
				self.base_steps[bp_idx[i]].tangent = np.asarray(final_tan[i])


	def calculate_angle_bw_tangents(self, base_step, masked=False):
		"""
		To calculate angle (Radian) between two tangent vectors of smoothed helical axis.

		**Arguments:**
				* ``base_step (1D list)``: List of two base-steps for which angle will be calculated.

					*Example:*
						
						``[5, 50]          # Calculate angle between tangent vectors of 5th and 50th base-steps``
			
				* ``masked (bool)``: ``Dfault=False``: To skip specific frames/snapshots. dnaMD.DNA.mask array should be set to use this functionality. This array contains boolean (either ``True`` or ``False``) value for each frame to mask the frames. Presently, mask array is automatically generated during :meth:`dnaMD.DNA.generate_smooth_axis` method to skip those frames where 3D fitting curve was not successfull within the given critera.

		**Returns:**
			* ``angle (1D array)``: Array of calculated angle of length is equal to number of frames. When ``masked`` is applied, length of this array can be smaller than total number of frames.

		"""

		if (len(base_step) != 2):
			raise ValueError("See, documentation for step usage!!!")
		
		tangent1, idx1 = self.get_parameters('Helical axis tangent', bp=[base_step[0]], bp_range=False, masked=masked)	
		tangent2, idx2 = self.get_parameters('Helical axis tangent', bp=[base_step[1]], bp_range=False, masked=masked)

		angle = []
		for i in range(len(tangent1[0])):
			tmp_angle = vector_angle(tangent1[0][i], tangent2[0][i])
			angle.append(tmp_angle)

		return np.asarray(angle)

def dev_bps_vs_parameter(dnaRef, bpRef, dnaSubj, bpSubj, parameter, err_type='std', bp_range=True, merge_bp=1, merge_method='mean'):
	"""To calculate deviation in the given parameters of a Subject DNA with respect to a Reference DNA along the base-pairs/steps.
	
		*Deviation = Reference_DNA(parameter) - Subject_DNA(parameter)*

		.. warning:: Number of base-pairs/steps should be similar in reference and subject DNA.
		
	**Arguments:**
	
		* ``dnaRef  (DNA object)``:   Reference DNA
		
		
		* ``bpRef (1D list) or (1D array)``: base-pairs or base-steps to consider from Reference DNA
				Example: ::
					
					bp = [6]                                # bp_range = False
					bp = [4,15]                             # bp_range = True
					bp = range(4,15)                        # bp_range = False
					bp = np.arange(4,15)                    # bp_range = False

		
		* ``dnaSubj (DNA object)``:   Subject DNA. Number of base-pairs in Reference and Subject DNA **should be** same.
		
		
		* ``bpSubj (1D list) or (1D array)``: base-pairs or base-steps to consider from Reference DNA. Foe more, see above example of ``bpSubj``.
		
		
		* ``parameter (string)``:   Name of a base-pair or base-step or helical base-step parameter
									For details about accepted keywords, see ``parameter`` in the method :meth:`DNA.get_parameters`.
		
		* ``error_type   (string)``:  Method of error estimation.
				Currently accepted method as follows:
				
				* ``error = 'std'``   : Standard Deviation
				* ``error = 'acf'``   : Standard error using autocorrelation time (requires: g_analyze)
				* ``error = 'block'`` : Standard error using block averaging method (requires: g_analyze)

			.. warning::
					to calculate errors by using ``error = 'acf'`` or ``error = 'block'`` , GROMACS tool ``g_analyze`` should be present in ``$PATH``.
			
		* ``bp_range (bool)``: Shown above. if True, bp should be range otherwise list or numpy array
	    
		* ``merge_bp  (int)``: Number of base-pairs or steps to merge for creating the small DNA segments
		
		* ``merge_method  (string)``: Method to calculate the parameter of a DNA segment from local parameters of all base-pairs/steps that are between the range given through ``bp``. 
				Currently accepted keywords are as follows:

				* ``merge_method = mean``: Average of local parameters
				* ``merge_method = sum``: Sum of local parameters
            
			

	**Returns:**
		* ``bpRef       (1D array)``: base-pair/step numbers of reference DNA. If ``merge_bp>1``, middle number will is returned.`
		* ``bpSubj      (1D array)``: base-pair/step numbers of subject DNA. If ``merge_bp>1``, middle number will is returned.`
		* ``deviation   (1D array)``: Deviation in the parameter of subject DNA with respect to reference DNA.
		* ``error       (1D array)``: Standard error of respective deviation
	"""

	bpRef, RefAvgValue, RefError = dnaRef.get_mean_error(bpRef, parameter, err_type = err_type, bp_range=True, merge_bp=merge_bp, merge_method=merge_method)

	bpSubj, SubjAvgValue, SubjError = dnaSubj.get_mean_error(bpSubj, parameter, err_type = err_type, bp_range=True, merge_bp=merge_bp, merge_method=merge_method)
	
	if len(bpRef) != len(bpSubj):
		raise ValueError("Number (%d) of bp/bps/segments in reference DNA does not match with the number (%d) of subject DNA." % (len(bpRef), len(bpSubj)))
		exit(1)

	deviation = RefAvgValue-SubjAvgValue
	error = np.sqrt((RefError**2)+(SubjError**2))

	return bpRef, bpSubj, deviation, error

def dev_parameters_vs_axis(dnaRef, dnaSubj, parameter, bp, axis ='Z', bp_range=True, windows=10, err_type='block'):
	"""To calculate deviation in the given parameters of a Subject DNA to Reference DNA along the given axis.
	
		*Deviation = Reference_DNA(parameter) - Subject_DNA(parameter)*
		
	**Arguments:**
	
		* ``dnaRef  (DNA object)``:   Reference DNA
		
		* ``dnaSubj (DNA object)``:   Subject DNA. Number of base-pairs in Reference and Subject DNA **should be** same.
		
		* ``parameter (string)``:   Name of a base-pair or base-step or helical base-step parameter
									For details about accepted keywords, see ``parameter`` in the method :meth:`DNA.get_parameters`.
		
		* ``bp (1D list) or (1D array)``: base-pairs or base-steps to analyze
				Example: ::
					
					bp = [6]                                # bp_range = False
					bp = [4,15]                             # bp_range = True
					bp = range(4,15)                        # bp_range = False
					bp = np.arange(4,15)                    # bp_range = False
					bp = [2,5,6,7,9,12,18]                  # bp_range = False


		* ``bp_range (bool)``: Shown above. if True, bp should be range otherwise list or numpy array
	    
		* ``windows  (int)``: Number of bins along the axis

	**Returns:**
		* ``deviation       (1D array)``: length = no. of windows; Deviation in the parameter for two given DNA
		* ``deviation_error (1D array)``: length = no. of windows; Standard error in deviation fo each window/bin
		* ``axis         (1D array)``: length = no. of windows; average position of window/bin along given axis
		* ``axis_error   (1D array)``: length = no. of windows; Standard error in average position of window/bin along given axis
	"""
	RefParam, ref_bp_idx = dnaRef.get_parameters(parameter,bp,bp_range)
	RefAxis, dummy = dnaRef.get_parameters('Helical {0}-axis' .format(axis),bp,bp_range)
	SubjParam, subj_bp_idx = dnaSubj.get_parameters(parameter,bp,bp_range)

	mean_axis = np.mean(RefAxis,axis=1)
	meanRefParam = np.mean(RefParam,axis=1)
	meanSubjParam = np.mean(SubjParam,axis=1)
	
	maxAxis = np.amax(mean_axis)
	minAxis = np.amin(mean_axis)
	axis_range = (maxAxis-minAxis)/windows

	Ref_param_error = get_error(dnaRef.time, RefParam, len(ref_bp_idx), err_type=err_type)
	Ref_axis_error = get_error(dnaRef.time, RefAxis, len(ref_bp_idx), err_type=err_type)
	subj_param_error = get_error(dnaSubj.time, SubjParam, len(subj_bp_idx), err_type=err_type)

	merged_ref_param =[]
	merged_subj_Param = []
	merged_Ref_param_error = []
	merged_Ref_axis_error = []
	merged_subj_param_error = []

	final_axis = []

	for i in range(windows):
		start = minAxis+(i*axis_range)
		end = start + axis_range
		idx = []
		for j in range(len(mean_axis)):
			if((start <= mean_axis[j]) and (end>mean_axis[j])):
				idx.append(j)
		if(len(idx)>0):
			merged_ref_param.append(meanRefParam[idx])
			merged_subj_Param.append(meanSubjParam[idx])
			final_axis.append( start + (end-start)/2 )
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
		
		final_ref_axis_error.append(np.sqrt( (merged_Ref_axis_error[i]**2).sum() ))
		final_ref_param_error.append(np.sqrt( (merged_Ref_param_error[i]**2).sum() ))
		final_subj_param_error.append(np.sqrt( (merged_subj_param_error[i]**2).sum() ))

	deviation,error = get_deviation(final_ref_param,final_ref_param_error,final_subj_param,final_subj_param_error)
	
	return deviation, error, final_axis, final_ref_axis_error

	

def get_error(time,x,sets,err_type='block'):
	"""To estimate error using block averaging method
	
	.. warning::
				It requires ``g_analyze`` of GROMACS package. ``g_analyze`` should be present in ``$PATH``.

	**Arguments:**
	
		* ``time (1D list) or (1D array)``: time
		* ``x	 (2D list) or (2D array)``: Shape of (nset, nframe); where *nset* is number of set and *nframe* is total number of frames. *nframe* should be equal to length of time list/array
		* ``sets (int)``: Number of sets (*nset*)
		* ``err_type (string)``: Error estimation by autocorrelation method ``err_type='acf'`` or block avearaging method ``err_type='block'``
	
	**Return:**
			``error (1D array)``: Of length = number of sets (*nset*)

	"""
	for i in range(sets):
		if (len(time) != len(x[i])):
			print '\nError: number of frame in time {0} mismatched with {1} of x[{2}]!!\n' .format(len(time),len(x[i]),i)
			exit(1)
	if not((err_type=='block') or (err_type=='acf')):
		print '\nWarning: Method {0} is not implemented. Swtiching to \'acf\'.\n' .format(err_type)
		err_type = 'acf'

	error = []
	char_set = string.ascii_lowercase
	name = ''.join(random.sample(string.ascii_lowercase, 10))
	
	filename = name+'.xvg'
	eefile = 'ee_'+name+'.xvg'
	acfile = 'acf_'+name+'.xvg'

	fout = open(filename,'w')
	for i in range(len(time)):
		fout.write('{0}' .format(time[i]))
		for j in range(sets):
			fout.write('\t{0}' .format(x[j][i]))
		fout.write("\n")
	fout.close()

	command = 'g_analyze -f {0} -ee {1} -ac {2} -fitfn exp' .format(filename,eefile,acfile)
	p = sub.Popen(command.split(), stdout=sub.PIPE, stderr=sub.PIPE)
	out, outputerror = p.communicate()
	lines = out.split('\n')

	if (err_type=='block'):
		for line in lines:
			if(re.match('Set',line)):
				temp = line.split()
				error.append(float(temp[3]))

	if (err_type=='acf'):
		acf_time = []
		for line in lines:
			if(re.match('COR: Correlation time',line)):
				temp = line.split('=')
				acf_time.append(abs(float(temp[1].split()[0])))
	
		total_time = float(time[-1])-float(time[0])
		dt = total_time/len(time)
		for i in range(sets):
			if(acf_time[i]>=dt):
				n_indp = total_time/acf_time[i]
				tmp_err = np.std(x[i])/np.sqrt(n_indp)
			else:
				tmp_err = np.std(x[i])/np.sqrt(len(time))
			error.append(tmp_err)
	
	os.remove(filename)
	os.remove(eefile)
	os.remove(acfile)
	if os.path.isfile('fitlog.log'):
		os.remove('fitlog.log')
	return np.array(error)

def get_deviation(Ref,RefErr,x,xerr):
	if (len(Ref) != len (x)):
		print "\nErrori: Number of base pairs/steps mismatch from reference to target!!\n"
		exit(1)
	Ref = np.array(Ref)
	RefErr = np.array(RefErr)
	x = np.array(x)
	xerr = np.array(xerr)

	covariance = np.cov(Ref,x)

	deviation = Ref-x
	error = np.sqrt( (RefErr*RefErr) + (xerr*xerr))

	return deviation, error
	
def frenet_serret(xyz):
	r''' Frenet-Serret Space Curve Invariants

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
		mag=np.sum(xyz**2,axis=1)**0.5
		imag=np.where(mag==0)
		mag[imag]=np.finfo(float).eps

		if n>1:
			return np.tile(mag,(n,1)).T
			
		return mag.reshape(len(mag),1) 

	xyz = np.asarray(xyz)
	n_pts = xyz.shape[0]
	
	if n_pts == 0:
		raise ValueError('xyz array cannot be empty')
    
	dxyz=np.gradient(xyz)[0]
	ddxyz=np.gradient(dxyz)[0]
	
	#Tangent
	T=np.divide(dxyz,magn(dxyz,3))
	
	#Derivative of Tangent
	dT=np.gradient(T)[0]    
	
	#Normal
	N = np.divide(dT,magn(dT,3))    
	
	#Binormal
	B = np.cross(T,N)    
    
	#Curvature
	k = magn(np.cross(dxyz,ddxyz),1)/(magn(dxyz,1)**3)    
    
	#Torsion 
	#(In matlab was t=dot(-B,N,2))
	t = np.sum(-B*N,axis=1)
	#return T,N,B,k,t,dxyz,ddxyz,dT   
    
	
	return T,N,B,k,t


def read_data_file(FileName,cols_equal=True):
	infile = open(FileName,'r')
	data = []
	len_data = 0
	i=1
	for line in infile:
		line = line.rstrip('\n')
		if not line.strip():
			continue
		if(re.match('#|@',line)==None):
			temp = map(float,line.split())
			if(cols_equal):
				if (i==1):
					len_data = len(temp)
				if (len(temp) != len_data):
					print 'WARNING: Number of column mis match at line {0} in {1}; skipping remaining part\n' .format(i,FileName)
					break
				data.append(temp)
			i = i+1
	data = np.array(data).T
	return data

def get_idx_of_bp_parameters(bp, parameters, bp_range):
	param_idx = []
	if(bp_range):
		bp_idx = np.arange(bp[0]-1,bp[1])
	else:
		bp_idx = np.subtract(bp,1)

	if(len(parameters)!=0):
		#param_idx = np.hstack((np.subtract(parameters,1),[parameters[-1]]))
		param_idx = np.subtract(parameters,1)
	
	return bp_idx, param_idx

def read_param_file(FileName,parameters, bp, bp_range, word=False):

	sys.stdout.write("\nReading file : %s\n" % FileName)
	sys.stdout.flush()
	
	def get_frame_data(block, parameters, bp_idx):
		block = np.array(block).T
		temp_data = (block[parameters,:])[:,bp_idx].copy()
		return temp_data

	def get_time(line):
		dummy, temp_time = line.split('=')
		return temp_time 

	infile = open(FileName,'r')
	data = []
	time = []
	frame_number = 0

	bp_idx, param_idx = get_idx_of_bp_parameters(bp,parameters,bp_range)

	block = []
	for line in infile:

		#Removing last new line charecter
		line = line.rstrip('\n')
		
		#Skipping blank/empty line
		if not line.strip():
			continue
		
		#Getting Time tag and time => Starting of new frame
		if(re.match('# Time',line)!=None):
			
			if((frame_number<100) and (frame_number%10==0)) or ((frame_number<1000) and (frame_number%100==0)) or ((frame_number<10000) and (frame_number%1000==0)) or ((frame_number<100000) and (frame_number%10000==0)) or ((frame_number<1000000) and (frame_number%100000==0)):
				sys.stdout.write("\rReading frame %d" % frame_number)
				sys.stdout.flush()
			
			frame_number += 1
			
			#if(frame_number==5000):
				#break
			
			#Getting time
			time.append(get_time(line))
			
			#Getting parameters/values for base-pairs
			if(len(block)>0):
				data.append(get_frame_data(block,param_idx,bp_idx))
			block = []
			continue

		#Skipping other lines starting with '#' tag
		if(re.match('#',line)!=None):
			continue
		
		if not word:
			block.append(map(float,line.split()))
		else:
			temp = []
			split_line = line.split()
			for word in split_line:
				if word != '---':
					temp.append(float(word))
				else:
					temp.append(None)
			block.append(temp)
	
	#For last frame
	data.append(get_frame_data(block,param_idx,bp_idx))
	block = []
	data_transpose = np.array(data).T

	sys.stdout.write("\nFinishid reading.... Total number of frame read =  %d\n" % frame_number)
	sys.stdout.flush()

	return data_transpose, time

def distance(x, y):
	x = np.asarray(x)
	y = np.asarray(y)
	return np.linalg.norm(x-y)

def vector_angle(x,y):
	dot = np.dot(x,y)
	cross = np.cross(x,y)
	cross_modulus = np.sqrt((cross*cross).sum())
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

		if count>4:
			mask = False
			sys.stdout.write('\n|frame:{0:>10}| WARNING: Fitting failed with \"smooth = {1}\"; Trying with \"smooth = {2}\".....\n' .format(nframe, smooth, smooth+100))
			sys.stdout.flush()
			
			smooth = smooth+100
			count = 1
			
			del orig_x
			del orig_y
			del orig_z
			
			orig_x = RawX.copy()
			orig_y = RawY.copy()
			orig_z = RawZ.copy()
		
		points = fill_point * len(orig_x)
		
		nest = -1
		tckp,u = splprep([orig_x, orig_y, orig_z], s=smooth, k=spline, nest=nest)
		
		xnew, ynew,znew = splev(np.linspace(0,1, points), tckp)
		
		new_axis = np.array([ xnew, ynew, znew ]).T

		angle = []
		dist = []
		del_idx = []
		last_idx = len(orig_x)-1
		
		for nbp in range(len(bp_idx)):
			start = nbp * fill_point
			end = start + fill_point
			xsmooth.append(xnew[start:end].mean())
			ysmooth.append(ynew[start:end].mean())
			zsmooth.append(znew[start:end].mean())
	
		for j in range(1, len(xsmooth)-1):
			prev = np.array([xsmooth[j-1], ysmooth[j-1], zsmooth[j-1]])
			curr = np.array([xsmooth[j], ysmooth[j], zsmooth[j]])
			nex = np.array([xsmooth[j+1], ysmooth[j+1], zsmooth[j+1]])
			angle.append( math.degrees(vector_angle((prev-curr),(curr-nex))) )

		for j in range(1, len(orig_x)-1):
			prev = np.array([orig_x[j-1], orig_y[j-1], orig_z[j-1]])
			curr = np.array([orig_x[j], orig_y[j], orig_z[j]])
			nex = np.array([orig_x[j+1], orig_y[j+1], orig_z[j+1]])
			dist.append(distance(prev, curr)+distance(curr, nex))
				
		for j in range(len(angle)):
					
			if angle[j] > cut_off_angle and not angle[j] > (180-cut_off_angle):
				
				del orig_x
				del orig_y
				del orig_z
				bsmooth = False
				
				max_idx = np.argsort(dist)[::-1]

				for k in range(count):
					del_idx.append(max_idx[k]+1)
					del_idx.append(max_idx[k]+2)
					if max_idx[k]==0:
						del_idx.append(0)
					if max_idx[k]==last_idx:
						del_idx.append(last_idx)

				del_idx = list(set(del_idx))
						
				sys.stdout.write('\r|frame:{0:>10}| WARNING: Bending angle [{1}-{2}-{3}] = {4:.2f} is more than cut-off angle {5};\n                     Four maximum distances between three adjacent axis positions = ({6[0]:.1f}, {6[1]:.1f}, {6[2]:.1f}, {6[3]:.1f});\n                     Deleting {7} original helical axis positions to remove possible fitting artifact...\n' .format(nframe, j-1, j, j+1, angle[j], cut_off_angle, np.sort(dist)[::-1], del_idx))
				sys.stdout.flush()
				mask = True

				if smooth>=10000:
					sys.stdout.write('\n\n|frame:{0:>10}| WARNING: Maximum Bending Angle = {1:.2f} at index {2}, which might be artifect. Please, check by visualizing PDB trajectory file...\n\n' .format(nframe, angle[j], j))
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
