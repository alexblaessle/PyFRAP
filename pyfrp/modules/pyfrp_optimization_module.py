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

"""Optimization module for PyFRAP toolbox.

Currently contains all functions necessary to transform a constrained FRAP optimization problem into
a unconstrained one, making it suitable to Nelder-Mead optimization algorithm. 

"""

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#Numpy/Scipy
import numpy as np

#PyFRAP
from . import pyfrp_fit_module 

from .pyfrp_term_module import *

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

def constrObjFunc(x,fit,debug,ax,returnFit):
	
	"""Objective function when using Constrained Nelder-Mead.
	
	Calls :py:func:`pyfrp.modules.pyfrp_optimization_module.xTransform` to transform x into
	constrained version, then uses :py:func:`pyfrp.modules.pyfrp_fit_module.FRAPObjFunc` to 
	find SSD.
	
	Args:
		x (list): Input vector, consiting of [D,(prod),(degr)].
		fit (pyfrp.subclasses.pyfrp_fit): Fit object.
		debug (bool): Display debugging output and plots.
		ax (matplotlib.axes): Axes to display plots in.
		returnFit (bool): Return fit instead of SSD.
	
	Returns:
		 float: SSD of fit. Except ``returnFit==True``, then will return fit itself. 
	"""
	
	
	LBs, UBs = buildBoundLists(fit)
	
	x=xTransform(x,LBs,UBs)

	ssd=pyfrp_fit_module.FRAPObjFunc(x,fit,debug,ax,returnFit)
	
	return ssd

def xTransform(x,LB,UB):
	
	"""Transforms ``x`` into constrained form, obeying upper 
	bounds ``UB`` and lower bounds ``LB``.
	
	.. note:: Will add tiny offset to LB(D), to avoid singularities.
	
	Idea taken from http://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd--fminsearchcon
	
	Args:
		x (list): Input vector, consiting of [D,(prod),(degr)].
		LB (list): List of lower bounds for ``D,prod,degr``.
		UB (list): List of upper bounds for ``D,prod,degr``.
	
	Returns:
		list: Transformed x-values. 
	"""
	
	#Make sure everything is float
	x=np.asarray(x,dtype=np.float64)
	LB=np.asarray(LB,dtype=np.float64)
	UB=np.asarray(UB,dtype=np.float64)
	
	#Check if LB_D==0, then add a little noise to it so we do not end up with xtrans[D]==0 and later have singularities when scaling tvec
	if LB[0]==0:
		LB[0]=1E-10
	
	#Determine number of parameters to be fitted
	nparams=len(x)

	#Make empty vector
	xtrans = np.zeros(np.shape(x))
	
	# k allows some variables to be fixed, thus dropped from the
	# optimization.
	k=0

	for i in range(nparams):

		#Upper bound only
		if UB[i]!=None and LB[i]==None:
		
			xtrans[i]=UB[i]-x[k]**2
			k=k+1
			
		#Lower bound only	
		elif UB[i]==None and LB[i]!=None:
			
			xtrans[i]=LB[i]+x[k]**2
			k=k+1
		
		#Both bounds
		elif UB[i]!=None and LB[i]!=None:
			
			xtrans[i] = (np.sin(x[k])+1.)/2.*(UB[i] - LB[i]) + LB[i]
			xtrans[i] = max([LB[i],min([UB[i],xtrans[i]])])
			k=k+1
		
		#No bounds
		elif UB[i]==None and LB[i]==None:
		
			xtrans[i] = x[k]
			k=k+1
			
		#Note: The original file has here another case for fixed variable, but since we made the decision earlier which when we call frap_fitting, we don't need this here.
	
	return xtrans	
		
def transformX0(x0,LB,UB):
	
	"""Transforms ``x0`` into constrained form, obeying upper 
	bounds ``UB`` and lower bounds ``LB``.
	
	Idea taken from http://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd--fminsearchcon
	
	Args:
		x0 (list): Input initial vector, consiting of [D,(prod),(degr)].
		LB (list): List of lower bounds for ``D,prod,degr``.
		UB (list): List of upper bounds for ``D,prod,degr``.
	
	Returns:
		list: Transformed x-values. 
	"""
	
	x0u = list(x0)
	
	nparams=len(x0)
	
	k=0
	for i in range(nparams):
		
		#Upper bound only
		if UB[i]!=None and LB[i]==None:
			if UB[i]<=x0[i]:
				x0u[k]=0
			else:
				x0u[k]=sqrt(UB[i]-x0[i])	
			k=k+1
			
		#Lower bound only
		elif UB[i]==None and LB[i]!=None:
			if LB[i]>=x0[i]:
				x0u[k]=0
			else:
				x0u[k]=np.sqrt(x0[i]-LB[i])	
			k=k+1
		
		
		#Both bounds
		elif UB[i]!=None and LB[i]!=None:
			if UB[i]<=x0[i]:
				x0u[k]=np.pi/2
			elif LB[i]>=x0[i]:
				x0u[k]=-np.pi/2
			else:
				x0u[k] = 2*(x0[i] - LB[i])/(UB[i]-LB[i]) - 1;
				#shift by 2*pi to avoid problems at zero in fminsearch otherwise, the initial simplex is vanishingly small
				x0u[k] = 2*np.pi+np.arcsin(max([-1,min(1,x0u[k])]));
			k=k+1
		
		#No bounds
		elif UB[i]==None and LB[i]==None:
			x0u[k] = x[i]
			k=k+1
	
	return x0u

def buildBoundLists(fit):
	
	"""Builds list of lower bounds and upper bounds.
	
	Args:
		fit (pyfrp.subclasses.pyfrp_fit): Fit object.
		
	Returns:
		tuple: Tuple containing: 
		
			* LBs (list): List of lower bounds.
			* UBs (list): List of upper bounds.
			
	
	
	"""
	
	LBs=[fit.LBD]+int(fit.fitProd)*[fit.LBProd]+int(fit.fitDegr)*[fit.LBDegr]+len(fit.ROIsFitted)*[fit.LBEqu]
	UBs=[fit.UBD]+int(fit.fitProd)*[fit.UBProd]+int(fit.fitDegr)*[fit.UBDegr]+len(fit.ROIsFitted)*[fit.UBEqu]
	
	return LBs,UBs