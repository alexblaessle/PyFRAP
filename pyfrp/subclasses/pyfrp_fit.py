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

"""Essential PyFRAP module containing :py:class:`pyfrp.subclasses.pyfrp_fit.fit` class. 
"""


#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

#Numpy/Scipy
import numpy as np

#PyFRAP Modules
from pyfrp.modules import pyfrp_misc_module 
from pyfrp.modules import pyfrp_plot_module
from pyfrp.modules import pyfrp_fit_module
from pyfrp.modules import pyfrp_stats_module


from pyfrp.modules.pyfrp_term_module import *

#Time 
import time

#===========================================================================================================================================================================
#Class definitions
#===========================================================================================================================================================================
	
class fit:
	
	"""Main fit class of PyFRAP. 
	
	The purpose of the fit class is to save all attributes used for fitting PyFRAP simulation results to 
	data analysis results. The main attributes are:
	
		* Fitting algorithm specifics:
		
			* Fitting algorithm, see also :py:func:`setOptMeth`.
			* Stopping criteria, see also :py:func:`setMaxfun` and :py:func:`setOptTol`.
			* Initial guess, see also :py:func:`getX0`.
			* Boundaries, see also :py:func:`getBounds`.
		
		* Fitting options:
		
			* ``fitProd``, see also :py:func:`getFitProd`.
			* ``fitDegr``, see also :py:func:`getFitDegr`.
			* ``fitPinned``, see also :py:func:`getFitPinned`.
			* ``equOn``, see also :py:func:`getEqu`.
			* ``fitCutOffT``, see also :py:func:`getFitCutOffT`.
			
		* The ROIs to be fitted, see also :py:func:`getROIsFitted`.
		
		* Fitting results, see also :py:func:`printResults`.
		
		* Fitted vectors.
		
	The most important methods are:
	
		* :py:func:`run`: Runs fitting.
		* :py:func:`addROIByName`: Adds ROI to be used for fitting.
		* :py:func:`getX0`: Builds and returns current initial guess.
		* :py:func:`getBounds`: Builds and returns current bounds.
		* :py:func:`computeStats`: Compares post-fitting statistics.
		
	The fit uses simulation and data vectors stored in all :py:class:`pyfrp.subclasses.pyfrp_ROI.ROI` objects defined in 
	``ROIsFitted`` list to compute the optimal values for ``DOptMu`` (``prodOpt`` or ``degrOpt`` if ``fitProd`` or ``fitDegr``
	is selected, respectively).  
			
	After calling :py:func:`run`, will automatically compute proper ``x0`` via :py:func:`getX0` and :py:func:`getBounds`.	
	
	Args:
		embryo (pyfrp.subclasses.pyfrp_embryo.embryo): Embryo object that fit belongs to.
		name (str): Name of fit.
	
	"""

	
	#Create new fit object
	def __init__(self,embryo,name):
		
		#General Settings
		self.name=name
		self.embryo=embryo
		
		#Optimization algorithm settings
		self.optMeth="Constrained Nelder-Mead"
		self.maxfun=1000
		self.optTol=1e-10
		
		#Dataseries selection
		self.ROIsFitted=[]
		
		#What parameters to fit
		self.fitProd=False
		self.fitDegr=False
		
		#Equalization and pinning
		self.equOn=True
		self.fitPinned=True
		self.equFacts=[]
		self.LBEqu=0.1
		self.UBEqu=3.
		
		#Intial guess
		self.x0=[10,0,0.]
		
		#Bounds
		self.LBProd=0.
		self.UBProd=100.
		self.LBDegr=0.
		self.UBDegr=100.
		self.LBD=0.01
		self.UBD=300.
		self.bounds=None
		
		#More settings
		self.kineticTimeScale=1.
		self.bruteInitD=False		
		
		#Cutting tvec option
		self.fitCutOffT=False
		self.cutOffT=150
		self.cutOffStepSim=self.embryo.simulation.stepsSim
		self.cutOffStepData=self.embryo.nFrames
		
		#Fit tracking
		self.saveTrack=0
		self.trackedParms=[]
		self.trackedFits=[]
		
		#Fitted Vectors
		self.fittedVecs=[]
		self.tvecFit=None
		self.dataVecsFitted=[]
		
		#Results
		self.SSD=10000000
		self.DOptMu=None
		self.DOptPx=None
		self.prodOpt=None
		self.degrOpt=None
		self.success=None
		self.iterations=None
		self.fcalls=None
		
		#Statistics
		self.Rsq=None
		self.MeanRsq=None
		self.RsqByROI={}
		
		#Empty result dataseries
		self.tvecFit=embryo.tvecData

	def addROI(self,r):
		
		"""Adds ROI to the list of fitted ROIs.
		
		Args:
			r (pyfrp.subclasses.pyfrp_ROI.ROI): ROI to be used for fitting.
		
		Returns:
			list: Updated list of ROIs used for fitting.
			
		"""
		
		if r not in self.ROIsFitted:
			self.ROIsFitted.append(r)
			self.x0.append(1.)
			
		return self.ROIsFitted
	
	def addROIByName(self,name):
		
		"""Adds ROI to the list of fitted ROIs, given a specific name.
		
		Args:
			name (str): Name of ROI to be used for fitting.
		
		Returns:
			list: Updated list of ROIs used for fitting.
			
		"""
		
		r=self.embryo.getROIByName(name)
		return self.addROI(r)
	
	def addROIById(self,Id):
		
		"""Adds ROI to the list of fitted ROIs, given a specific ROI Id.
		
		Args:
			Id (int): Id of ROI to be used for fitting.
		
		Returns:
			list: Updated list of ROIs used for fitting.
			
		"""
		
		r=self.embryo.getROIById(Id)
		return self.addROI(r)
	
	def getROIsFitted(self):
		
		"""Returns list of ROIs used for fitting.
		
		Returns:
			list: list of ROIs used for fitting.
		
		"""
		
		return self.ROIsFitted
		
	
	def removeROI(self,r):
		
		"""Removes ROI from the list of fitted ROIs.
		
		Args:
			r (pyfrp.subclasses.pyfrp_ROI.ROI): ROI to be removed.
		
		Returns:
			list: Updated list of ROIs used for fitting.
			
		"""
		
		if r in self.ROIsFitted:
			
			idx=self.ROIsFitted.index(r)
			self.x0.pop(3+idx)
			self.ROIsFitted.remove(r)
			
		return self.ROIsFitted
	
	def getX0(self):
		
		"""Returns initial guess of fit in the form that is useful for 
		the call of the optimization algorithm.
		
		Copies x0 into local variable to pass to solver, pop entries that are 
		currently not needed since they are turned off via ``fitProd`` or ``fitDegr``.
		
		Always appends initial guess for equalization factors, even though they might not been used.
		
		.. note:: Always gets executed at the start of ``run``.
		
		Returns:
			list: Currently used x0.
			
		"""
		
		x0=list(self.x0)	
		if self.fitProd and self.fitDegr:
			pass	
		elif self.fitProd and  not self.fitDegr:	
			x0.pop(2)
		elif not self.fitProd and self.fitDegr:
			x0.pop(1)
		elif not self.fitProd and not self.fitDegr:
			x0=list(self.x0)
			x0.pop(2)
			x0.pop(1)
		
		return x0
	
	def reset2DefaultX0(self):
		
		"""Resets initial guess x0 to its default form.
		
		The default form of x0 is 
		
		>>> [10., 0. ,0. , 1.,1.,1.]
		
		The last entries are the initial guess of equlalization factors and is
		set to be list of of ones of the same length of ``ROIsFitted``.
		
		Returns:
			list: New initial guess x0.
		
		"""
		
		equFacts=len(self.ROIsFitted)*[1.]
		
		self.x0=[10,0,0]+equFacts
		
		return self.x0
	
	def getBounds(self):
		
		"""Generates tuple of boundary tuples, limiting parameters 
		varied during SSD minimization.
		
		Will generate exactly the boundary tuple that is currently 
		useful to the optimization algorithm, meaning that only 
		values that are needed  since they are turned on via ``fitProd`` or ``fitDegr``
		will be included into tuple.
		
		Will use values that are stored in ``LBx`` and ``UBx``, where ``x`` is 
		``D``, ``Prod``, or ``Degr`` for the creation of the tuples. 
		
		Will also add a tuple of bounds defined via ``LBEqu`` and ``UBEqu`` for each 
		ROI in ``ROIsFitted``.
		
		.. note:: Always gets executed at the start of ``run``.
		
		Returns:
			tuple: Boundary value tuple.
			
		"""
	
		if self.fitProd and self.fitDegr:
			bnds = [(self.LBD, self.UBD), (self.LBD, self.UBD),(self.LBD,self.UBD)]
			ranges=[slice(self.LBD,self.UBD,1),slice(self.LBProd,self.UBProd,10),slice(self.LBDegr,self.UBDegr,10)]
		elif self.fitProd and  not self.fitDegr:	
			bnds = [(self.LBD, self.UBD), (self.LBProd, self.UBProd)]
			ranges=[slice(self.LBD,self.UBD,1),slice(self.LBProd,self.UBProd,10)]
		elif not self.fitProd and self.fitDegr:
			bnds = [(self.LBD, self.UBD), (self.LBDegr, self.UBDegr)]
			ranges=[slice(self.LBD,self.UBD,1),slice(self.LBDegr,self.UBDegr,10)]
		elif not self.fitProd and not self.fitDegr:
			bnds = [(self.LBD, self.UBD),]
			ranges=[1,self.UBD]
		
		bnds=bnds+len(self.ROIsFitted)*[(self.LBEqu,self.UBEqu)]
		ranges=ranges+len(self.ROIsFitted)*[slice(self.LBEqu,self.UBEqu,0.2)]
		
		if self.optMeth=='brute':
			self.bounds=tuple(ranges)
		else:
			self.bounds=tuple(bnds)
			
		return self.bounds
	
	def resultsToVec(self):
		
		"""Puts results back in vector as optimization algorithm would return it.
		
		Returns:
			list: Result vector.
		"""
		
		x=[self.DOptPx]
		
		if self.fitProd:
			x.append(self.prodOpt)	
		if self.fitDegr:	
			x.append(self.degrOpt)
		if self.equOn:
			x=x+list(self.equFacts)
		
		return x
	
	def getFittedParameterNames(self):
		
		"""Returns names of parameters that are selected for fitting.
		
		Returns:
			list: Names of parameters fitted.
			
		"""
		
		x=["DOptPx"]
		
		if self.fitProd:
			x.append("prod")	
		if self.fitDegr:	
			x.append("degr")
		if self.equOn:
			for r in self.ROIsFitted:
				x.append(r.name+" equFact")
			
		return x
	
	def run(self,debug=False,ax=None):
		
		"""Runs fit.
		
		Fitting is done by passing fit object to :py:func:`pyfrp.modules.pyfrp_fit_module.FRAPFitting`.
		This function then calls all necessary methods of fit to prepare it for optimization and
		then passes it to optimization algorithm.
		
		.. note:: If ``bruteInitD`` is turned on, will execute :py:func:`runBruteInit` instead.
		
		Keyword Args:
			debug (bool): Print debugging messages.
			ax (matplotlib.axes): Axes to show debugging plots in.
		
		Returns:
			pyfrp.subclasses.pyfrp_fit.fit: ``self``.
		
		"""
		
		if self.bruteInitD:
			self.runBruteInit(debug=debug,ax=ax)
		else:	
			self=pyfrp_fit_module.FRAPFitting(self,debug=debug,ax=ax)
		
		return self
	
	def runBruteInit(self,debug=False,ax=None,steps=5,x0Ds=[]):
		
		"""Runs fit for different initial guesses of the diffusion constant D, then
		selects the one that actually yielded the minimal SSD.
		
		Initially guesses are generated with :py:func:`getBruteInitDArray` if no array ``x0Ds``
		is given.
		
		Fitting is done by passing fit object to :py:func:`pyfrp.modules.pyfrp_fit_module.FRAPFitting`.
		This function then calls all necessary methods of fit to prepare it for optimization and
		then passes it to optimization algorithm.
		
		Will select the initial guess that yielded the minimal SSD and then rerun with this x0 again, making
		sure that everything is updated in fit object.
		
		Keyword Args:
			debug (bool): Print debugging messages.
			ax (matplotlib.axes): Axes to show debugging plots in.
			steps (int): How many initial guesses to generate.
			x0Ds (list): Array with possible initial guesses for D.
		
		Returns:
			pyfrp.subclasses.pyfrp_fit.fit: ``self``.
		
		"""
		
		if x0Ds==[]:
			x0Ds=self.getBruteInitDArray(steps=steps)
		SSDs=[]
		
		for x0D in x0Ds:
			
			if debug:
				print("Trying x0(D) = ", x0D)
			
			self.setX0D(x0D)
			self=pyfrp_fit_module.FRAPFitting(self,debug=debug,ax=ax)
			
			SSDs.append(self.SSD)
			
		idxOpt=SSDs.index(min(SSDs))
		
		if debug:
			print("x0(D) yielding best result = ", x0Ds[idxOpt]) 
		
		self.setX0D(x0Ds[idxOpt])
		self=pyfrp_fit_module.FRAPFitting(self,debug=debug,ax=ax)
		
		return self
		
	def getBruteInitDArray(self,steps=5):
		
		"""Generates array of different possibilities to be used as initial guess 
		for D.
		
		If ``LBD`` and ``UBD`` is given, will simply divide the range between the two in 4 equidistant values.
		Otherwise will vary around ``x0`` in 2 orders of magnitude.
		
		Keyword Args:
			steps (int): How many initial guesses to generate.
		
		Returns:
			list: Array with possible initial guesses for D.
		
		"""
		
		if self.LBD!=None:
			LB=self.LBD
		else:
			LB=0.01*self.getX0D()
			
		if self.UBD!=None:
			UB=self.UBD
		else:
			UB=100*self.getX0D()
		
		x0=np.linspace(LB+1E-10,UB-1E-10,steps)
		
		return list(x0)
	
	def assignOptParms(self,res):
		
		r"""Assigns optimal parameters found by optimization algorithm to
		attributes in fit object depending on fit options chosen.
		
		Args:
			res (list): Result array from optimization algorithm.
		
		Returns:
			tuple: Tuple containing:
			
				* DOptPx (float): Optimal diffusion coefficient in :math:`\frac{\mathrm{px}^2}{s}}`.
				* prod (float): Optimal production rate in :math:`\frac{\[c\]}{s}}`.
				* degr (float): Optimal degradation rate in :math:`\frac{1}{s}}`.
				* DOptMu (float): Optimal diffusion coefficient in :math:`\frac{\mu\mathrm{m}^2}{s}}`.
				
		"""
		
		if self.fitProd and self.fitDegr:
			self.DOptPx=res[0]
			self.prodOpt=res[1]/self.kineticTimeScale
			self.degrOpt=res[2]/self.kineticTimeScale
			
		elif self.fitProd and not self.fitDegr:
			self.DOptPx=res[0]
			self.prodOpt=res[1]/self.kineticTimeScale
			self.degrOpt=self.x0[2]/self.kineticTimeScale
			
		elif not self.fitProd and self.fitDegr:
			self.DOptPx=res[0]
			self.prodOpt=self.x0[1]/self.kineticTimeScale
			self.degrOpt=res[1]/self.kineticTimeScale
			
		elif not self.fitProd and not self.fitDegr:
			self.DOptPx=res[0]
			self.prodOpt=self.x0[1]/self.kineticTimeScale
			self.degrOpt=self.x0[2]/self.kineticTimeScale
		
		
		self.DOptMu=self.DOptPx*self.embryo.convFact**2
		
		return self.DOptPx, self.prodOpt, self.degrOpt, self.DOptMu
	
	def plotFit(self,ax=None,legend=True,title=None,show=True):
		
		"""Plots fit, showing the result for all fitted ROIs.
		
		.. note:: If no ``ax`` is given, will create new one.
		
		.. image:: ../imgs/pyfrp_fit/fit.png
		
		Keyword Args:
			ax (matplotlib.axes): Axes used for plotting.
			legend (bool): Show legend.
			title (str): Title of plot.
			show (bool): Show plot.
		
		Returns:
			matplotlib.axes: Axes used for plotting.
		
		"""
		
		for r in self.ROIsFitted:
			ax=r.plotFit(self,ax=ax,legend=legend,title=title,show=show)
			
		return ax
	
	def printResults(self):
		
		"""Prints out main results of fit."""
		
		printObjAttr('name',self)
		printObjAttr('DOptMu',self)
		printObjAttr('DOptPx',self)
		printObjAttr('prodOpt',self)
		printObjAttr('degrOpt',self)
		printObjAttr('equFacts',self)
		printObjAttr('success',self)
		printObjAttr('Rsq',self)
		printObjAttr('MeanRsq',self)
		printObjAttr('RsqByROI',self)
		
		return True
	
	def printAllAttr(self):
		
		"""Prints out all attributes of fit object.""" 
		
		printAllObjAttr(self)
	
	def resultsToDict(self):
		
		"""Extracts all important results into dictionary, making
		it easier for printout or csv extraction.
		
		"""
		
		parms=["DOptMu","DOptPx","prodOpt","degrOpt","success","Rsq","MeanRsq","fitDegr","fitProd","fitPinned","equOn","x0"]
		
		dic=pyfrp_misc_module.objAttr2Dict(self,attr=parms)
		
		roiNames=pyfrp_misc_module.objAttrToList(self.ROIsFitted,"name")
		#equFacts=np.asarray(self.equFacts).astype(str)
		
		dic["ROIsFitted"]=" , ".join(roiNames)
		
		
		for i in range(len(roiNames)):
			
			if self.equOn:
				dic["equFactor "+roiNames[i]]=self.equFacts[i]
			else:
				dic["equFactor "+roiNames[i]]=""
			
			dic["Rsq("+roiNames[i]+")"]=self.RsqByROI[roiNames[i]]
				
		return dic
		
	
	def setBruteInitD(self,b):
		
		"""Turns on/off if the initial guess of for the diffusion rate D should be bruteforced.
		
		Args:
			b (bool): Flag value.
			
		Returns:
			bool: Current flag value.
		
		"""
		
		self.bruteInitD=b
		
		return self.bruteInitD
	
	def setOptMeth(self,m):
		
		"""Sets optimization method.
		
		Available optimization methods are:
		
			* Constrained Nelder-Mead
			* Nelder-Mead
			* TNC
			* L-BFGS-B
			* SLSQP
			* brute
			* BFGS
			* CG
		
		See also http://docs.scipy.org/doc/scipy-0.17.0/reference/generated/scipy.optimize.minimize.html and
		http://docs.scipy.org/doc/scipy-0.17.0/reference/generated/scipy.optimize.brute.html#scipy.optimize.brute .
		
		You can find out more about the constrained Nelder-Mead algorithm in the documentation of 
		:py:func:`pyfrp.modules.pyfrp_optimization_module.constrObjFunc`.
		
		Args:
			m (str): New method.
			
		"""
		
		self.optMeth=m
		return self.optMeth
	
	def getOptMeth(self):
		
		"""Returns the currently used optimization algorithm.
		
		Returns:
			str: Optimization algorithm.
			
		"""
		
		return self.optMeth
		
	def isFitted(self):
		
		"""Checks if fit already has been run and succeeded.
		
		Returns:
			bool: ``True`` if success.
			
		"""
		
		return self.DOptMu!=None
		
	def setEqu(self,b):
		
		"""Turns on/off equalization.
		
		Args:
			b (bool): New flag value.
			
		Returns:
			bool: New flag value.
		
		"""
		
		self.equOn=b
		return self.equOn
	
	def setFitPinned(self,b):
		
		"""Turns on/off if pinned series are supposed to be fitted.
		
		Args:
			b (bool): New flag value.
			
		Returns:
			bool: New flag value.
		
		"""
		
		self.fitPinned=b
		return self.fitPinned
	
	def setFitProd(self,b):
		
		"""Turns on/off if production is supposed to be considered in fit.
		
		Args:
			b (bool): New flag value.
			
		Returns:
			bool: New flag value.
		
		"""
		
		self.fitProd=b
		return self.fitProd
	
	def setFitDegr(self,b):
		
		"""Turns on/off if degradation is supposed to be considered in fit.
		
		Args:
			b (bool): New flag value.
			
		Returns:
			bool: New flag value.
		
		"""
		
		self.fitDegr=b
		return self.fitDegr
	
	def setSaveTrack(self,b):
		
		"""Turns on/off if fitting process is supposed to be stored.
		
		This then can then be used to following the convergence of 
		the optimization algorithm and possibly to identify local minima.
		
		Args:
			b (bool): New flag value.
			
		Returns:
			bool: New flag value.
		
		"""
		
		self.saveTrack=b
		return self.saveTrack
	
	def setFitCutOffT(self,b):
		
		"""Turns on/off if only a certain fraction of the timeseries
		is supposed to be fitted.
		
		.. warning:: This option is currently VERY experimental. Fitting might
		   crash.
		
		Args:
			b (bool): New flag value.
			
		Returns:
			bool: New flag value.
		
		"""
		
		
		printWarning("CutOffT Option is currently VERY experimental. Fitting might crash.")
		self.fitCutOffT=b
		return self.fitCutOffT
	
	def setCutOffT(self,t):
		
		self.cutOffT=t
		return self.cutOffT
	
	def setMaxfun(self,m):
		
		"""Sets maximum number of function evaluations at 
		which optimization algorithm stops.
		
		Args:
			m (int): New maximum number of function evaluations.
		
		"""
		
		self.maxfun=m
		return self.maxfun
	
	def setOptTol(self,m):
		
		"""Sets tolerance level at which optimization algorithm stops.
		
		Args:
			m (float): New tolerance level.
		
		"""
		
		self.optTol=m
		return self.optTol
	
	def setLBD(self,b):
		
		"""Sets the lower bound for the diffusion rate.
		
		Args:
			b (float): New lower bound for diffusion rate.
		
		"""
		
		self.LBD=b
		return self.LBD
	
	def setLBProd(self,b):
		
		"""Sets the lower bound for the production rate.
		
		Args:
			b (float): New lower bound for production rate.
		
		"""
		
		self.LBProd=b
		return self.LBProd
	
	def setLBDegr(self,b):
		
		"""Sets the lower bound for the degradation rate.
		
		Args:
			b (float): New lower bound for degradation rate.
		
		"""
		
		self.LBDegr=b
		return self.LBDegr
	
	def setUBD(self,b):
		
		"""Sets the upper bound for the diffusion rate.
		
		Args:
			b (float): New upper bound for diffusion rate.
		
		"""
		
		self.UBD=b
		return self.UBD
	
	def setUBProd(self,b):
		
		"""Sets the upper bound for the production rate.
		
		Args:
			b (float): New upper bound for production rate.
		
		"""
		
		self.UBProd=b
		return self.UBProd
	
	def setUBDegr(self,b):
		
		"""Sets the upper bound for the degradation rate.
		
		Args:
			b (float): New upper bound for degradation rate.
		
		"""
		
		self.UBDegr=b
		return self.UBDegr
			
	def getEqu(self):
		
		"""Returns equalization flag.
			
		Returns:
			bool: Current flag value.
		
		"""
		
		return self.equOn
	
	def getFitPinned(self):
		
		"""Returns flag controlling if pinned timeseries are supposed to be used 
		for fitting.
			
		Returns:
			bool: Current flag value.
		
		"""
		
		return self.fitPinned
	
	def getFitProd(self):
		
		"""Returns flag controlling if a production term is supposed to be used 
		for fitting.
			
		Returns:
			bool: Current flag value.
		
		"""
		
		return self.fitProd
	
	def getFitDegr(self):
		
		"""Returns flag controlling if a degredation term is supposed to be used 
		for fitting.
			
		Returns:
			bool: Current flag value.
		
		"""
		
		
		return self.fitDegr
	
	def getSaveTrack(self):
			
		"""Returns flag controlling if whole fitting process is supposed to be saved
		in fit object.
			
		Returns:
			bool: Current flag value.
		
		"""
		
		return self.saveTrack
	
	def getFitCutOffT(self):
			
		"""Returns flag controlling if only the first ``cutOffT`` timesteps are supposed to be fitted.
		
		.. warning:: This option is currently VERY experimental. Fitting might
		   crash.
		
		Returns:
			bool: Current flag value.
		
		"""
		
		return self.fitCutOffT
	
	def getCutOffT(self):
		
		"""Returns timepoint at which timeseries are cut if ``fitCutOffT`` is turned on.
		
		.. warning:: This option is currently VERY experimental. Fitting might
		   crash.
		
		Returns:
			float: Timepoint.
		
		"""
		
		return self.cutOffT
	
	def getMaxfun(self):
		
		"""Returns maximum number of function evaluations at 
		which optimization algorithm stops.
		
		Returns:
			int: Current maximum number of function evaluations.
		
		"""
		
		return self.maxfun
	
	def getOptTol(self):
		
		"""Returns tolerance level at which optimization algorithm stops.
		
		Returns:
			float: Current tolerance level.
		
		"""
		
		return self.optTol
	
	def getLBD(self):
		
		"""Returns the lower bound for the diffusion rate.
		
		Returns:
			float: Current lower bound for diffusion rate.
		
		"""
		
		return self.LBD
	
	def getLBProd(self):
		
		"""Returns the lower bound for the production rate.
		
		Returns:
			float: Current lower bound for production rate.
		
		"""
		
		
		return self.LBProd
	
	def getLBDegr(self):
		
		"""Returns the lower bound for the degradation rate.
		
		Returns:
			float: Current lower bound for degradation rate.
		
		"""
		
		return self.LBDegr
	
	def getUBD(self):
		
		"""Returns the upper bound for the diffusion rate.
		
		Returns:
			float: Current upper bound for diffusion rate.
		
		"""
		
		return self.UBD
	
	def getUBProd(self):
		
		"""Returns the upper bound for the production rate.
		
		Returns:
			float: Current upper bound for production rate.
		
		"""
		
		return self.UBProd
	
	def getUBDegr(self):
		
		"""Returns the upper bound for the degradation rate.
		
		Returns:
			float: Current upper bound for degradation rate.
		
		"""
		
		return self.UBDegr
	
	def setKineticTimeScale(self,s):
		
		"""Sets the kinetic time scale factor used for fitting.
		
		Args:
			s (float): New kinetic time scale factor.
			
		"""
		
		self.kineticTimeScale=s
		return self.kineticTimeScale
	
	def getKineticTimeScale(self):
		
		"""Returns the kinetic time scale factor used for fitting.
		
		Returns:
			float: Current kinetic time scale factor.
			
		"""
		
		return self.kineticTimeScale
	
	def setName(self,s):
			
		"""Sets name of fit.
		
		Args:
			s (str): New name of fit.
			
		"""
		
		self.name=s
		return self.name
	
	def getName(self):
		
		"""Returns name of fit.
		
		Returns:
			str: Name of fit.
			
		"""
		
		return self.name
	
	def setX0Equ(self,x):
		
		"""Sets the initial guess for the equalization factor.
		
		.. note:: Does this for all ROIs in ROIsFitted.
		
		Args:
			x (float): Initial guess for equalization factor.
		"""
		
		for i in range(3,len(self.x0)):
			self.x0[i]=x
		
		return self.x0
	
	def getX0Equ(self,x):
		
		"""Returns the initial guess for the equalization factor for
		all ROIs fitted.
		
		Returns:
			list: Initial guess for equalization factor.
		"""
		
		return self.x0[3:]
	
	def setX0D(self,x):
		
		"""Sets the initial guess for the diffusion rate.
		
		Args:
			x (float): Initial guess for diffusion rate.
		"""
		
		self.x0[0]=x
		return self.x0[0]
	
	def setX0Prod(self,x):
		
		"""Sets the initial guess for the production rate.
		
		Args:
			x (float): Initial guess for production rate.
		"""
		
		self.x0[1]=x
		return self.x0[1]
	
	def setX0Degr(self,x):
		
		"""Sets the initial guess for the degradation rate.
		
		Args:
			x (float): Initial guess for degradation rate.
		"""
		
		self.x0[2]=x
		return self.x0[2]
	
	def getX0D(self):
		
		"""Returns the initial guess for the diffusion rate.
		
		Returns:
			float: Initial guess for diffusion rate.
		"""
		
		return self.x0[0]
	
	def getX0Prod(self):
		
		"""Returns the initial guess for the production rate.
		
		Returns:
			float: Initial guess for production rate.
		"""
		
		return self.x0[1]
	
	def getX0Degr(self):
		
		"""Returns the initial guess for the degradation rate.
		
		Returns:
			float: Initial guess for degration rate.
		"""
		
		return self.x0[2]
	
	def setX0(self,x):
		
		"""Sets the initial guess ``x0``.
		
		Argument ``x`` needs to have length 3, otherwise it is being rejected. 
		
		.. note:: If ``fitProd`` or ``fitDegr`` are not chosen, the values in
		   ``x0`` are going to be used as static parameters.
		
		
		Args:
			x (list): New desired initial guess.
			
		Returns:
			list: New initial guess.
			
		
		"""
		
		if len(x)==3:
			self.x0=x
		else:
			printError("Length of x0 is not 3, not going to change it.")
				
		return self.x0
	
	def checkPinned(self):
		
		"""Checks if all ROIs in ``ROIsFitted`` have been pinned.
		
		Returns:
			bool: ``True`` if all ROIs have been pinned, ``False`` else.
		"""
		
		b=True
		for i,r in enumerate(self.ROIsFitted):
			b = b + len(self.embryo.tvecData)==len(r.dataVecPinned) + len(self.embryo.simulation.tvecSim)==len(r.simVecPinned) 
		return b
	
	def checkSimulated(self):
		
		"""Checks if all ROIs in ``ROIsFitted`` have been simulated.
		
		Returns:
			bool: ``True`` if all ROIs have been simulated, ``False`` else.
		"""
		
		b=True
		for r in self.ROIsFitted:
			b = b + len(r.simVec)==len(self.embryo.simulation.tvecSim)
		return b
		
	def updateVersion(self):
		
		"""Updates fit object to current version, making sure that it possesses
		all attributes.
		
		Creates a new fit object and compares ``self`` with the new fit object.
		If the new fit object has a attribute that ``self`` does not have, will
		add attribute with default value from the new fit object.
		
		
		Returns:
			pyfrp.subclasses.pyfrp_fit.fit: ``self``
			
		"""
		
		fittemp=fit(self.embryo,"temp")
		pyfrp_misc_module.updateObj(fittemp,self)
		return self
	
	def computeStats(self):
		
		"""Computes stastics for fit.
		
		Statistics include:
			
			* ``MeanRsq``
			* ``Rsq``
			* ``RsqByROI``
		
		"""
		
		self=pyfrp_stats_module.computeFitRsq(self)
		
	def printRsqByROI(self):
		
		"""Prints out Rsq value per ROI.
		"""
		
		print("Rsq Values by ROI for fit ", self.name)
		printDict(self.RsqByROI)
		
	def getNParmsFitted(self,inclEqu=True):
		
		"""Returns the number of parameters fitted in this fit.
		
		.. note:: If equlalization is turned on, each ROI in ``ROIsFitted``
		   counts as an extra parameter.
		   
		Example: We fit production and equalization for 2 ROIs, then we have fitted
		
			* D
			* degradation
			* equalization ROI 1
			* equalization ROI 2
			
		leading to in total 4 fitted parameters.
			
		Keyword Args:
			inclEqu (bool): Include equalization as additional fitted parameter.
		
		Returns:
			int: Number of parameters fitted.
		
		"""
		
		return 1+int(self.getFitProd())+int(self.getFitDegr())+int(inclEqu)*int(self.getEqu())*len(self.equFacts)
	
	def plotLikehoodProfiles(self,epsPerc=0.1,steps=100,debug=False):
		
		"""Plots likelihood profiles for all fitted parameters.
		
		.. warning:: Since we don't yet fit the loglikelihood function, we only plot the 
		   SSD. Even though the SSD is proportional to the loglikelihood, it should be used
		   carefully.
		   
		See also :py:func:`pyfrp.modules.pyfrp_fit_module.plotFitLikehoodProfiles`.   
		
		Keyword Args:
			epsPerc (float): Percentage of variation.
			steps (int): Number of values around optimal parameter value.
			debug (bool): Show debugging messages
		
		Returns:
			list: List of matplotlib.axes objects used for plotting.
		
		"""
		
		axes=pyfrp_fit_module.plotFitLikehoodProfiles(self,epsPerc=0.1,steps=100,debug=debug)
		
		return axes
	
	def setOpts(self,opts):
		
		"""Sets a list of options.
		
		Options are given as a dictionary and then subsequentially set.
		
		Args:
			opts (dict): Options.
			
			
		
		"""
		
		for opt in opts:
			try:
				setattr(self,opt,opts[opt])
			except AttributeError:
				printError("Cannot set fit option " + opt +". Option does not exist.")
				
				
		
		