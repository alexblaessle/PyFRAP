#=====================================================================================================================================
#Copyright
#=====================================================================================================================================

#Copyright (C) 2014 Alexander Blaessle, Patrick Mueller and the Friedrich Miescher Laboratory of the Max Planck Society
#This software is distributed under the terms of the GNU General Public License.

#This file is part of PyFRAP.

#PyFRAP is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

#===========================================================================================================================================================================
#Module Description
#===========================================================================================================================================================================

"""Integration module for PyFRAP toolbox. 

.. warning:: Might get merged with simulation module at some point.

"""

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================


from numpy import linalg as LA

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

def getAvgConc(val,cvs,ind):
	
	"""Integrates simulation result over specific set of indices.
	
	Args:
		val (fipy.CellVariable): PDE solution variable.
		cvs (numpy.ndarray): Array containing cell volumes.
		ind (list): List of indices.
	
	Returns:
		float: Integration result.
	
	"""
	
	if len(ind)>0:
		if hasattr(val,'value'):
			return sum(val.value[ind]*cvs[ind])/sum(cvs[ind])
		else:
			return sum(val[ind]*cvs[ind])/sum(cvs[ind])
	else:
		return 0.

def calcTetSidelengths(point0,point1,point2,point3):

	"""Calculates sidelengths of tetrahedron given by 4 points.
	
	.. note:: Taking point0 as base point.
	
	"""

	vec1=point1-point0
	vec2=point2-point0
	vec3=point3-point0
	
	norms=[LA.norm(vec1),LA.norm(vec2), LA.norm(vec3)]
	
	return norms
