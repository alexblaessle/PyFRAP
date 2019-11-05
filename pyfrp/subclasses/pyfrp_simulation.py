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

#Simulation class for PyFRAP toolbox, including following classes:

#(1) simulation

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#Numpy/Scipy
import numpy as np
import scipy.interpolate as interp 

#PyFRAP classes
from . import pyfrp_mesh

#PyFRAP modules
from pyfrp.modules import pyfrp_plot_module
from pyfrp.modules import pyfrp_sim_module
from pyfrp.modules import pyfrp_img_module
from pyfrp.modules import pyfrp_idx_module
from pyfrp.modules import pyfrp_misc_module
from pyfrp.modules.pyfrp_term_module import *

#Plotting
import matplotlib.pyplot as plt

#itertools 
import itertools

#===========================================================================================================================================================================
#Class definitions
#===========================================================================================================================================================================

class simulation(object):
	
	"""PyFRAP simulation class. 
	
	Stores all important properties about how FRAP simulation is performed, such as:

	"""
	
	
	def __init__(self,embryo):
		
		#Naming/ID
		self.embryo=embryo
		
		#IC image
		self.ICimg=None
		
		self.restoreDefaults()
		
	def restoreDefaults(self):
		
		"""Restores default parameters for simulations."""
		
		#PDE specific
		self.D=50.
		self.prod=0.0
		self.degr=0.0
		
		#Mesh specific
		self.mesh=self.setMesh(pyfrp_mesh.mesh(self))
		
		#IC specific
		self.ICmode=3
		self.IC=None
		
		#Only necessary if ICmode=4 (ideal ICs)
		self.bleachedROI=None
		self.valOut=None
		
		#Time specific
		self.stepsSim=3000
		self.tvecSim=np.linspace(self.embryo.tStart,self.embryo.tEnd,self.stepsSim)
		
		#Save simulation
		self.saveSim=False
		self.vals=[]
		
		#Solver details
		self.solver="PCG"
		self.iterations=1000
		self.tolerance=1E-10
	
	def setSolver(self,solver):
		
		"""Sets solver to use.
		
		Implemented solvers are:
		
			* PCG
			* LU
			
		Args:
			solver (str): Solver to use.
		
		Returns:
			str: Current solver
		
		"""
		
		if solver not in ["PCG","LU"]:
			printWarning("Unknown solver " + solver +". This might lead to problems later")
		
		self.solver=solver
		return self.solver
	
	def getSolver(self):
		
		"""Returns current solver.
		
		Returns:
			str: Current solver
		
		"""
	
		return self.solver
	
	def setTolerance(self,tol):
		
		"""Sets tolerance of solver.
			
		Args:
			tol (float): New tolerance.
		
		Returns:
			float: Current tolerance.
		
		"""
	
		self.tolerance=tolerance
		return self.tolerance
	
	def getTolerance(self):
		
		"""Returns current tolerance.
		
		Returns:
			str: Current solver
		
		"""
		
		
		return self.tolerance
	
	def setIterations(self,tol):
		
		"""Sets iterations of solver.
			
		Args:
			tol (float): New iterations.
		
		Returns:
			float: Current iterations.
		
		"""
	
		self.iterations=iterations
		return self.iterations
	
	def getIterations(self):
		
		"""Returns current iterations.
		
		Returns:
			str: Current solver
		
		"""
		
		return self.iterations
	
	def setICMode(self,m):
		
		"""Sets the mode of initial conditions.
		
		Initial condition modes are defined as:
		
			* 0: **ROI-based**: Mesh nodes get assigned the value of the first entry ``dataVec`` of 
			  the ROI covering them. Note: If a mesh node is covered by two ROIs, will assign the value
			  of the ROI that is last in embryo's ``ROIs`` list. See also 
			  :py:func:`pyfrp.modules.pyfrp_sim_module.applyROIBasedICs`.
			* 1: **Radial**: Concentrations get radially approximated around some center.
			  See also :py:func:`pyfrp.modules.pyfrp_sim_module.applyRadialICs`.
			* 2: **Imperfect**: Interpolates ``ICimg`` onto mesh, fills nodes that are not covered
			  by image with ``embryo.analysis.concRim`` and then mimics imperfect bleaching via sigmoid
			  function in z-direction.
			  See also :py:func:`pyfrp.modules.pyfrp_sim_module.applyImperfectICs`.
			* 3: **Inpterpolated**: Interpolates ``ICimg`` onto mesh, fills nodes that are not covered
			  by image with ``embryo.analysis.concRim``.
			  See also :py:func:`pyfrp.modules.pyfrp_sim_module.applyInterpolatedICs`.
			* 4: **Ideal**: Ideal ICs, that is, single value inside and single value outside of bleached region.
			  See also :py:func:`pyfrp.modules.pyfrp_sim_module.applyIdealICs`, :py:func:`setValOut` 
			  and :py:func:`setBleachedROI`.
		
		.. note:: The default mode is **Interpolated** (``m=3``) and is highly recommended to obtain most realistic results.
		
		Args:
			m (int): Which mode to be used.
			
		Returns:
			int: Current initial condition mode used.
		
		"""
		
		if m not in list(range(5)):
			printError("ICmode = " +  m + " is not defined. Not going to change ICmode" )
			return self.ICmode
			
		self.ICmode=m
		
		if m==4:
			if self.bleachedROI==None:
				printWarning("bleachedROI is not set yet. This might lead to problems later.")
			if self.valOut==None:
				printWarning("valOut is not set yet. This might lead to problems later.")
		
		return self.ICmode
	
	def setBleachedROI(self,r):
		
		"""Sets bleached ROI that is used when ideal ICs (ICmode=4) is selected.
		
		Args:
			r (pyfrp.subclasses.pyfrp_ROI.ROI): ROI to be set bleached ROI.
			
		"""
		
		self.bleachedROI=r
		
		return self.bleachedROI
	
	def setValOut(self,v):
		
		"""Sets valOut that is used when ideal ICs (ICmode=4) is selected.
		
		Args:
			v (float): Value that is to assigned outside of bleachedROI.
		
		"""
		
		self.valOut=v
		
		return self.valOut
	
	def run(self,signal=None,embCount=None,showProgress=True,debug=False):
		
		"""Runs simulation.
		
		Checks if ROI indices are computed, if not, computes them. Then passes simulation
		object to :py:func:`pyfrp.modules.pyfrp_sim_module.simulateReactDiff`.
		
		Keyword Args:
			signal (PyQt5.QtCore.pyqtSignal): PyQT signal to send progress to GUI.
			embCount (int): Counter of counter process if multiple datasets are simulated. 
			debug (bool): Print debugging messages and show debugging plots.
			showProgress (bool): Print out progress.
		
		Returns:
			bool: True if success, False otherwise.
			
		
		"""
		
		self.updateVersion()
		
		if not self.embryo.checkROIIdxs()[1]:
			self.embryo.computeROIIdxs()
			
		if not isinstance(self.ICimg,np.ndarray) and self.ICmode in [2,3]:
			printWarning("No ICimg was specified, but it is required for selected ICmode="+str(self.ICmode)+". Will grab first image in "+self.embryo.fnDatafolder)
			try:
				self.setICimgByFn(self.embryo.getDataFolder()+'/'+self.embryo.getFileList()[0])
			except:
				printError("Was not able to set new ICimg. Will abort.")
				return False
			
		pyfrp_sim_module.simulateReactDiff(self,signal=signal,embCount=embCount,showProgress=showProgress,debug=debug)
		return True
	
	def setMesh(self,m):
		
		"""Sets mesh to a new mesh object.
	
		Args:
			m (pyfrp.subclasses.pyfrp_mesh.mesh): PyFRAP mesh object.
			
		Returns:
			pyfrp.subclasses.pyfrp_mesh.mesh: Updated mesh instance.
			
		
		"""
		
		self.mesh=m
		return self.mesh
	
	def setICimgByFn(self,fn):
		
		"""Sets image for initial condition interpolation given a filepath.
		
		Args:
			fn (str): Path to file.
			
		Returns:
			numpy.ndarray: New ICimg.
		
		"""
		
		img=pyfrp_img_module.loadImg(pyfrp_misc_module.fixPath(fn),self.embryo.dataEnc)
		
		return self.setICimg(img)
		
	def setICimg(self,img):
		
		"""Sets image for initial condition interpolation.
		
		Args:
			img (numpy.ndarray): A 2D image.
			
		Returns:
			numpy.ndarray: New ICimg.
		
		"""
		
		self.ICimg=img
		return self.ICimg
	
	def getICimg(self):
			
		"""Returns image for initial condition interpolation.
		
		Returns:
			numpy.ndarray: Current ICimg.
		
		"""
		
		
		return self.ICimg
	
	def showIC(self,ax=None,roi=None,nlevels=25,vmin=None,vmax=None,typ='contour',scale=True):
		
		"""Plots initial conditions applied to mesh in 2D or 3D.
		
		If ``roi`` is given, will only plot initial conditions for nodes inside ROI, else 
		will plot initial condition for all nodes in mesh.
		
		.. note:: Simulation needs to be run first before this plotting function
		   can be used.
		   
		Example:
		
		>>> simulation.plotIC(typ='contour')
		
		will produce the following:
		
		.. image:: ../imgs/pyfrp_simulation/showIC.png
		
		See also :py:func:`pyfrp.modules.pyfrp_plot_module.plotSolutionVariable` and :py:func:`pyfrp.subclasses.pyfrp_ROI.plotSolutionVariable`.
		
		Keyword Args:
			roi (pyfrp.subclasses.pyfrp_ROI.ROI): A PyFRAP ROI object.
			vmin (float): Overall minimum value to be displayed in plot.
			vmax (float): Overall maximum value to be displayed in plot.
			ax (matplotlib.axes): Axes used for plotting.
			nlevels (int): Number of contour levels to display.
			typ (str): Typ of plot.
			scale (bool): Equal axis in case of contour plot.
		
		Returns:
			matplotlib.axes: Axes used for plotting.
		
		"""
		
		
		if self.IC!=None:
			
			if ax==None:
				if typ=='surface':
					fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["IC"],proj=['3d'],sup="",tight=False)
				else:	
					fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["IC"],sup="",tight=False)
				ax=axes[0]
			
			if roi==None:
		
				x,y,z=self.getCellCenters()
			
				if vmin==None:
					vmin=min(self.IC)
				if vmax==None:
					vmax=max(self.IC)
				
				levels=np.linspace(vmin,1.01*vmax,nlevels)
				
				if typ=='contour':
					ax.tricontourf(x,y,self.IC,vmin=vmin,vmax=vmax,levels=levels)
					if scale:
						ax.autoscale(enable=True, axis='both', tight=True)
				elif typ=='surface':
					ax.plot_trisurf(x,y,self.IC,cmap='jet',vmin=vmin,vmax=vmax)
				else:
					printError("Unknown plot type "+ typ)
				
				ax.get_figure().canvas.draw()
		
			else:
				ax=roi.plotSolutionVariable(self.IC,ax=ax,nlevels=nlevels,vmin=vmin,vmax=vmax,typ=typ)
			
			return ax
			
		else:
			printWarning("IC is not generated yet. Run simulation first.")
			return None	
		
	def showICimg(self,ax=None,typ='contour',colorbar=True,scale=True,nlevels=25,vmin=None,vmax=None):
		
		"""Plots image used for initial condition either as contour or surface plot.
		
		.. image:: ../imgs/pyfrp_simulation/showICimg.png
		
		Keyword Args:
			ax (matplotlib.axes): Axes used for plotting.
			scale (bool): Equal axis.
			vmin (float): Overall minimum value to be displayed in plot.
			vmax (float): Overall maximum value to be displayed in plot.
			nlevels (int): Number of contour levels to display.
			
		Returns:
			matplotlib.axes: Axes used for plotting.
		
		"""
		
		#Check of entered plot type makes sense
		if typ not in ['contour','surface']:
			printError("Unknown plot type "+ typ)
			return ax
		
		if self.ICimg!=None:
			
			if vmin==None:
				vmin=min(self.ICimg.flatten())
			if vmax==None:
				vmax=max(self.ICimg.flatten())
			
			levels=np.linspace(vmin,1.01*vmax,nlevels)
			
			if ax==None:
				if typ=='surface':
					fig,axes = pyfrp_plot_module.makeSubplot([1,1],proj=['3d'],sup="",tight=False)
				else:
					fig,axes = pyfrp_plot_module.makeSubplot([1,1],sup="",tight=False)
				ax=axes[0]
				
			res=self.ICimg.shape[0]
			if 'quad' in list(self.embryo.analysis.process.keys()):
				X,Y=np.meshgrid(np.arange(res,2*res),np.arange(res,2*res))
			else:
				X,Y=np.meshgrid(np.arange(res),np.arange(res))
			
			if typ=='contour':
				plt_ICs=ax.contourf(X,Y,self.ICimg,levels=levels,vmin=vmin,vmax=vmax)
				if scale:
					plt.axis('equal')
			elif typ=='surface':
				plt_ICs=ax.plot_surface(X,Y,self.ICimg,cmap='jet',vmin=vmin,vmax=vmax)
				
			if colorbar:
				cb=plt.colorbar(plt_ICs,orientation='horizontal',pad=0.05,shrink=0.9)
			
			plt.draw()
			
			return ax
			
		else:
			printWarning("ICimg is not analyzed yet. Run data analysis first.")
			return None
	
	def computeInterpolatedICImg(self):
		
		"""Computes interpolation of initial condition image.
		
		Interpolation is done as in :py:func:`pyfrp.modules.pyfrp_sim_module.applyInterpolatedICs`.
		
		Returns:
			tuple: Tuple containing:
			
				* xInt (numpy.ndarray): Meshgrid x-coordinates.
				* yInt (numpy.ndarray): Meshgrid y-coordinates.
				* f (numpy.ndarray): Interpolated image.
	
		"""
		
		#Get image resolution and center of geometry
		res=self.ICimg.shape[0]
		center=self.embryo.geometry.getCenter()
		
		#Define x/y coordinates of interpolation
		if 'quad' in list(self.embryo.analysis.process.keys()):
			#Shift everything by center to fit with the mesh
			xInt = np.arange(center[0]+1, center[0]+res+1, 1)
			yInt = np.arange(center[1]+1, center[1]+res+1, 1)		
		else:
			xInt = np.arange(1, res+1, 1)
			yInt = np.arange(1, res+1, 1)
		
		#Generate interpolation function
		f=interp.RectBivariateSpline(xInt, yInt, self.ICimg, bbox=[None, None, None, None], kx=3, ky=3, s=0)
		
		return xInt, yInt, f
	
	def computeInterpolatedSolutionToImg(self,vals,roi=None,method='linear',res=None):
		
		"""Interpolates solution back onto 2D image.
		
		Uses ``scipy.interpolate.griddata``, see also http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.griddata.html
		
		If ``roi`` is specified, will only interpolate nodes of this ROI. 
		
		For more details about interpolation methods, check out 
		https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.griddata.html .
		
		Keyword Args:
			vals (numpy.ndarray): Solution to be interpolated.
			roi (pyfrp.subclasses.pyfrp_ROI.ROI): A PyFRAP ROI.
			method (str): Interpolation method.
			fillVal (float): Value applied outside of ROI.
			res (int): Resolution of resulting images in pixels.
			
		Returns:
			tuple: Tuple containing:
			
				* X (numpy.ndarray): Meshgrid x-coordinates.
				* Y (numpy.ndarray): Meshgrid y-coordinates.
				* interpIC (numpy.ndarray): Interpolated solution.
		
		"""
		
		#Get image resolution and center of geometry
		if res==None:
			res=self.ICimg.shape[0]
		
		#Build Empty Img
		X,Y=np.meshgrid(np.arange(res),np.arange(res))
		
		#Get cellcenters
		x,y,z=self.getCellCenters()
		##print x
		##print y
		##print z
		##print roi.meshIdx
		#print max(roi.meshIdx)
		#print len(x)
		
		if roi!=None:
			xInt=x[roi.meshIdx]
			yInt=y[roi.meshIdx]
			val=vals[roi.meshIdx]
		else:
			xInt=x
			yInt=y
			val=vals
		
		interpIC=interp.griddata((xInt,yInt),val,(X,Y),method=method)
		
		return X,Y,interpIC
	
	def computeInterpolatedICImg(self,roi=None):
		
		"""Interpolates ICs back onto 2D image.
		
		Uses ``scipy.interpolate.griddata``, see also http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.griddata.html
		
		If ``roi`` is specified, will only interpolate nodes of this ROI. 
		
		Keyword Args:
			roi (pyfrp.subclasses.pyfrp_ROI.ROI): A PyFRAP ROI.
			
		Returns:
			tuple: Tuple containing:
			
				* X (numpy.ndarray): Meshgrid x-coordinates.
				* Y (numpy.ndarray): Meshgrid y-coordinates.
				* interpIC (numpy.ndarray): Interpolated ICs.
		
		"""
		
		X,Y,interpIC=self.computeInterpolatedSolutionToImg(self.IC,roi=roi)
		
		return X,Y,interpIC
		
	def showInterpolatedIC(self,ax=None,roi=None):
		
		"""Shows ICs interpolated back onto 2D image.
		
		If ``roi`` is specified, will only interpolate nodes of this ROI. 
		
		See also :py:func:`computeInterpolatedIC`.
		
		.. image:: ../imgs/pyfrp_simulation/showInterpolatedIC.png
		
		Keyword Args:
			roi (pyfrp.subclasses.pyfrp_ROI.ROI): A PyFRAP ROI.
			ax (matplotlib.axes): Axes to be used for plotting.
			
		Returns:
			matplotlib.axes: Axes used for plotting.
		
		"""
		
		X,Y,interpIC=self.computeInterpolatedIC(roi=roi)
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["Interpolated IC"],sup="simulation")
			ax=axes[0]
			
		ax.imshow(interpIC)
		
		return ax
	
	def showInterpolatedICImg(self,ax=None):
		
		"""Shows interpolation of initial condition image.
		
		See also :py:func:`computeInterpolatedICImg`.
		
		.. image:: ../imgs/pyfrp_simulation/showInterpolatedICimg.png
		
		Keyword Args:
			ax (matplotlib.axes): Axes to be used for plotting.
			
		Returns:
			matplotlib.axes: Axes used for plotting.
		
		"""
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["Interpolated Image"],sup="simulation")
			ax=axes[0]
			
		xInt, yInt, f=self.computeInterpolatedICImg()	
		
		#print np.shape(xInt)
		
		#raw_input()
		X,Y=np.meshgrid(xInt,yInt)
		imgInt=np.zeros(np.shape(X))
		
		for i in range(np.shape(X)[0]):
			for j in range(np.shape(Y)[0]):
				imgInt[i,j]=f(X[i,j],Y[i,j])
				
		ax.imshow(imgInt)	
		
		return ax
	
	def compareICInterpolation(self,axes=None,roi=None):
		
		"""Shows initial image, its interpolation, the resulting initial 
		condition and its interpolation back onto an image.
		
		See also :py:func:`showICimg`, :py:func:`showInterpolatedICImg`, 
		:py:func:`showIC`, :py:func:`showInterpolatedIC`.
		
		Will create new axes if necessary.
		
		.. warning:: Some images might be flipped due to plotting functions. Will be fixed in future version.
		
		.. image:: ../imgs/pyfrp_simulation/ICcompare.png
		
		Keyword Args:
			roi (pyfrp.subclasses.pyfrp_ROI.ROI): A PyFRAP ROI.
			axes (matplotlib.axes): List of axes of length 4.
			
		Returns:
			list: List of axes.
		
		"""
	
		if axes==None:
			fig,axes = pyfrp_plot_module.makeSubplot([2,2],titles=["Original Image","Interpolated Image","IC","Reinterpolated IC"],sup="simulation")
		
		self.showICimg(ax=axes[0])
		self.showInterpolatedICImg(ax=axes[1])
		self.showIC(ax=axes[2],roi=roi)
		self.showInterpolatedIC(ax=axes[3],roi=roi)
		
		for ax in axes:
			pyfrp_plot_module.redraw(ax)
		
		return axes	
	
	def setTEnd(self,T):
		
		"""Updates timevector of simulation to end at new time point.
		
		.. note:: Keeps scaling.
		
		Args:
			T (float): New end timepoint.
			
		Returns:
			numpy.ndarray: New simulation time vector.
		
		"""
		
		wasLog=self.isLogTimeScale()
			
		self.tvecSim=np.linspace(self.embryo.tStart,T,self.stepsSim)
		
		if wasLog:
			self.toLogTimeScale()
		
		return self.tvecSim
	
	def getOptTvecSim(self,maxDExpectedPx):
		
		r"""Generates time vector that is optimal to fit 
		experiments with expected diffusion coefficients up
		to ``maxDExpectedPx``.
		
		Basically computes how long a simulation needs to run in
		seconds to capture the dynamics of an experiment with diffusion
		coefficient of ``maxDExpectedPx``. Does this by setting end time point to
		
		.. math:: t_{\mathrm{end,sim}} = \frac{D_{\mathrm{max. exp.}}}{D_{\mathrm{sim}}} t_{\mathrm{end,data}}
		
		.. note:: Keeps time scaling.
		
		Args:
			maxDExpectedPx (float): Maximum expected diffusion coefficient.
			
		Returns:
			numpy.ndarray: New simulation time vector.
		
		"""
		
		wasLog=self.isLogTimeScale()
			
		self.tvecSim=np.linspace(self.embryo.tStart,maxDExpectedPx/self.D*self.embryo.tEnd,self.stepsSim)
		
		if wasLog:
			self.toLogTimeScale()
		
		return self.tvecSim
		
	def toLogTimeScale(self,spacer=1E-10):
		
		"""Converts time vector for simulation to logarithmic scale.
		
		Keyword Args:
			spacer (float): Small offset to avoid log(0).
			
		Returns:
			numpy.ndarray: New simulation time vector.
		
		"""
		
		self.tvecSim=self.tvecSim[0]+np.logspace(np.log10(spacer+self.tvecSim[0]), np.log10(self.tvecSim[-1]), self.stepsSim)-spacer
		return self.tvecSim
		
	def toLinearTimeScale(self):
		
		"""Converts time vector for simulation to linear scale.
	
		Returns:
			numpy.ndarray: New simulation time vector.
		
		"""
		
		self.tvecSim=np.linspace(self.tvecSim[0],self.tvecSim[-1],self.stepsSim)	
		return self.tvecSim
	
	def toDefaultTvec(self):
		
		"""Sets time vector for simulation to default range.
		
		Default range is given by ``tStart`` and ``tEnd`` in ``embryo`` object
		and is linearly scaled.
		
		Returns:
			numpy.ndarray: New simulation time vector.
		
		"""
		
		self.tvecSim=np.linspace(self.embryo.tStart,self.embryo.tEnd,self.stepsSim)
		return self.tvecSim
	
	def updateTvec(self):
		
		"""Updates time vector for simulation to match 
		experiment start and end time.
		
		Does not change scaling of time vector.
		
		Returns:
			numpy.ndarray: New simulation time vector.
		
		"""
		
		if (self.tvecSim[1]-self.tvecSim[0])==(self.tvecSim[-1]-self.tvecSim[-2]):
			self.toLinearTimeScale()
		else:
			self.toLogTimeScale()
		return self.tvecSim
	
	def setTimesteps(self,n):
		
		"""Sets number of simulation time steps and updates
		time vector.
		
		Args:
			n (int): New number of time steps.
			
		Returns:
			int: New number of time steps.
			
		"""
		
		self.stepsSim=int(n)
		self.updateTvec()
		return self.stepsSim
	
	def setD(self,D):
		
		"""Sets diffusion coefficient used for simulation.
		
		Args:
			D (float): New diffusion coefficient in :math:`\mu\mathrm{m}^2/s`.
			
		Returns:
			float: New diffusion coefficient in :math:`\mathrm{px}^2/s`.
			
		"""
		
		self.D=D
		return self.D
	
	def getD(self):
		
		"""Returns current diffusion coefficient used for simulation.
		
		Returns:
			float: Current diffusion coefficient in :math:`\mathrm{px}^2/s`.
			
		"""
		
		return self.D
	
	def setProd(self,prod):
		
		"""Sets production rate used for simulation.
		
		Args:
			prod (float): New production rate in :math:`1/s`.
			
		Returns:
			float: New production rate in :math:`1/s`.
			
		"""
		
		self.prod=prod
		return self.prod
	
	def getProd(self):
		
		"""Returns production rate used for simulation.
		
		Returns:
			float: Current production rate in :math:`1/s`.
			
		"""
		
		return self.prod
	
	def setDegr(self,degr):
		
		"""Sets degradation rate used for simulation.
		
		Args:
			prod (float): New degradation rate in :math:`1/[c]s`.
			
		Returns:
			float: New degradation rate in :math:`1/[c]s`.
			
		"""
		
		self.degr=degr
		return self.degr
	
	def getDegr(self):
		
		"""Returns degradation rate used for simulation.
			
		Returns:
			float: Current degradation rate in :math:`1/[c]s`.
			
		"""
		
		return self.degr
	
	def setSaveSim(self,b):
		
		"""Sets flag if simulation should be saved.
		
		Args:
			b (bool): New flag value.
		
		Returns:
			bool: Updated flag value.
		
		"""
		
		self.saveSim=b
		
		return self.saveSim
	
	def getSaveSim(self):
		
		"""Returns flag if simulation should be saved.
		
		Returns:
			bool: Current flag value.
		
		"""
		
		return self.saveSim
	
	def isLogTimeScale(self):
		
		"""Returns if time spacing of simulation is logarithmic.
		
		Returns:
			bool: Time spacing is logarithmic.
		
		"""
		
		#Note: We round here to 5 decimals since diff seems to make numerical mistakes.
		
		return not round(np.diff(self.tvecSim)[0],5)==round(np.diff(self.tvecSim)[-1],5)
	
	def plotICStack(self,ROIs,withGeometry=True,vmin=None,vmax=None,ax=None,colorbar=False):
		
		"""Plots a stack of the initial conditions in a given list of ROIs.
		
		Will automatically compute the direction in which ROI lies in the 3D space and
		reduce the ROI into this plane for contour plot.
		
		If ``vmin=None`` or ``vmax=None``, will compute overall maximum and minimum values
		over all ROIs.
		
		.. image:: ../imgs/pyfrp_simulation/ICstack.png
		
		
		Args:
			phi (fipy.CellVariable): Simulation solution variable (or numpy array).
			ROIs (list): List of :py:class:`pyfrp.subclasses.pyfrp_ROI.ROI` objects.
			
		Keyword Args:
			withGeometry (bool): Show geometry inside plot.
			vmin (float): Overall minimum value to be displayed in plot.
			vmax (float): Overall maximum value to be displayed in plot.
			ax (matplotlib.axes): Axes used for plotting.
			colorbar (bool): Display color bar.
		
		Returns:
			matplotlib.axes: Axes used for plotting.
		
		"""
		
		ax = self.plotSolStack(self.IC,ROIs,withGeometry=withGeometry,vmin=vmin,vmax=vmax,ax=ax,colorbar=colorbar)
		
		return ax 
	
	
		
	def plotSolStack(self,phi,ROIs,withGeometry=True,vmin=None,vmax=None,ax=None,colorbar=False):
			
		"""Plots a stack of the solution variable in a given list of ROIs.
		
		Will automatically compute the direction in which ROI lies in the 3D space and
		reduce the ROI into this plane for contour plot.
		
		If ``vmin=None`` or ``vmax=None``, will compute overall maximum and minimum values
		over all ROIs.
		
		Args:
			phi (fipy.CellVariable): Simulation solution variable (or numpy array).
			ROIs (list): List of :py:class:`pyfrp.subclasses.pyfrp_ROI.ROI` objects.
			
		Keyword Args:
			withGeometry (bool): Show geometry inside plot.
			vmin (float): Overall minimum value to be displayed in plot.
			vmax (float): Overall maximum value to be displayed in plot.
			ax (matplotlib.axes): Axes used for plotting.
			colorbar (bool): Display color bar.
		
		Returns:
			matplotlib.axes: Axes used for plotting.
		
		"""
			
			
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["Simulation IC stack"],proj=['3d'])
			ax=axes[0]
		
		#Plot geometry
		if withGeometry:
			#self.embryo.geometry.updateGeoFile()
			ax=self.embryo.geometry.plotGeometry(ax=ax)
		
		#Find vmin/vmax over all ROIs
		vminNew=[]
		vmaxNew=[]
		for r in ROIs:
			vminNew.append(min(phi[r.meshIdx]))
			vmaxNew.append(max(phi[r.meshIdx]))
		
		if vmin==None:
			vmin=min(vminNew)	
		if vmax==None:
			vmax=min(vmaxNew)
			
		for r in ROIs:
			plane=r.getMaxExtendPlane()
			zs=r.getPlaneMidCoordinate()
			zdir=r.getOrthogonal2Plane()
			
			ax=r.plotSolutionVariable(phi,ax=ax,vmin=vmin,vmax=vmax,plane=plane,zs=zs,zdir=zdir,colorbar=colorbar)
				
		return ax
	
	def getSolutionVariableSmoothness(self,vals,roi=None):
		
		r"""Returns smoothness of solution variable. 
	
		Smoothness :math:`s` is computed as:
		
		.. math:: s=\frac{d_{\mathrm{max}}}{\bar{d}}
		
		where :math:`d_{\mathrm{max}}` is the maximum derivative from the nearest neighbour over the whole array, and 
		:math:`\bar{d}` the average derivation. Derivative from nearest neighbour is computed by
		
		.. math:: d=\frac{c-c_\mathrm{nearest}}{||\textbf{x}-\textbf{x}_\mathrm{nearest}||_2}
		
		.. note:: If ``roi!=None``, will only evaluate smoothness over this specific ROI.
		
		.. warning:: Nearest neighbour finding algorithm is slow. Should be changed to ``ckdTree`` at some point.
		
		Keyword Args:
			roi (pyfrp.subclasses.pyfrp_ROI.ROI): PyFRAP ROI.
		
		Returns:
			tuple: Tuple containing:
			
				* s (float): Smoothmess coefficient.
				* dmax(float): Maximum diff.
			
		"""
		
		#Retrieve values
		x,y,z=self.getCellCenters()
	
		if roi!=None:
			idxs=roi.meshIdx
		else:	
			idxs=list(range(len(x)))
		
		x=np.asarray(x)[idxs]
		y=np.asarray(y)[idxs]
		z=np.asarray(z)[idxs]
		vals=np.asarray(self.IC)[idxs]
		
		#Nearest neighbour
		idx,dist=pyfrp_idx_module.nearestNeighbour3D(x,y,z,x,y,z,k=2,minD=0)
		#Compute deriv
		diffs=[]
		for i in range(len(idx)):
			diffs.append(abs(vals[i]-vals[idx[i]]))
	
		#Get derivative
		deriv=np.asarray(diffs)/np.asarray(dist)	
				
		#Normalize
		s=max(deriv)/(np.mean(deriv))
		
		return s,max(deriv)
		
	def getICSmoothness(self,roi=None):
		
		"""Returns smoothness of initial condition. 
	
		See also :py:func:`getSolutionVariableSmoothness`.
		
		.. note:: If ``roi!=None``, will only evaluate smoothness over this specific ROI.
		
		Keyword Args:
			roi (pyfrp.subclasses.pyfrp_ROI.ROI): PyFRAP ROI.
		
		Returns:
			tuple: Tuple containing:
			
				* s (float): Smoothmess coefficient.
				* dmax(float): Maximum diff.
			
		"""
		
		return self.getSolutionVariableSmoothness(self.IC,roi=roi)
		
	def getICImgSmoothness(self):
		
		r"""Returns smoothness of initial condition. 
	
		See also :py:func:`pyfrp.modules.pyfrp_img_module.getICImgSmoothness`.
		
		Returns:
			tuple: Tuple containing:
			
				* s (float): Smoothmess coefficient.
				* dmax(float): Maximum diff.
		"""	
		
		return pyfrp_img_module.getImgSmoothness(self.ICimg)
		
	def rerun(self,signal=None,embCount=None,showProgress=True,debug=False):
		
		"""Reruns simulation.
		
		.. note:: Only works if simulation has been run before with ``saveSim`` enabled.
		
		See also :py:func:`pyfrp.modules.pyfrp_sim_module.rerunReactDiff`.
		
		Keyword Args:
			signal (PyQt5.QtCore.pyqtSignal): PyQT signal to send progress to GUI.
			embCount (int): Counter of counter process if multiple datasets are simulated. 
			debug (bool): Print debugging messages and show debugging plots.
			showProgress (bool): Print out progress.
		
		Returns:
			pyfrp.subclasses.pyfrp_simulation.simulation: Updated simulation instance.
			
		
		"""
		
		if not self.embryo.checkROIIdxs()[1]:
			self.embryo.computeROIIdxs()
		
		pyfrp_sim_module.rerunReactDiff(self,signal=signal,embCount=embCount,showProgress=showProgress,debug=debug)
		return True
		
	def updateVersion(self):
		
		"""Updates simulation object to current version, making sure that it possesses
		all attributes.
		
		Creates a new simulation object and compares ``self`` with the new simulation object.
		If the new simulation object has a attribute that ``self`` does not have, will
		add attribute with default value from the new simulation object.
		
		Returns:
			pyfrp.subclasses.pyfrp_simulation.simulation: ``self``
			
		"""
		
		simtemp=simulation(self.embryo)
		pyfrp_misc_module.updateObj(simtemp,self)
		return self
		
	def printAllAttr(self):
		
		"""Prints out all attributes of embryo object.""" 
		
		print("Simulation of embryo ", self.embryo.name, " Details.")
		printAllObjAttr(self)
		
	def mapOntoImgs(self,tvec=None,roi=None,fnOut="",showProgress=True,method='linear',fillVal=0.,scale=True,enc='uint16',res=None):
		
		"""Maps simulation solution back onto images.
		
		See also :py:func:`computeInterpolatedSolutionToImg`.
		
		.. note:: Only works if simulation has been run before and saved via ``saveSim``.
		
		For more details about interpolation methods, check out 
		https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.griddata.html .
		
		Keyword Args:
			tvec (numpy.ndarray): Timepoints at which solution is saved to image.
			roi (pyfrp.subclasses.pyfrp_ROI.ROI): PyFRAP ROI.
			fnOut (str): Path where images should be saved.
			showProgress (bool): Show progress of output.
			method (str): Interpolation method.
			fillVal (float): Value applied outside of ROI.
			scale (bool): Scale to encoding range.
			enc (str): Encoding.
			res (int): Resolution of resulting images in pixels.
		
		Returns:
			bool: True if everything ran through, False else.
		"""
		
		#Check if anything has been saved
		if len(self.vals)==0:
			printError("Simulation hasn't been saved, will not do anything.")
			return False
			
		#Grab tvec	
		if tvec==None:
			tvec=self.tvecSim
		
		#Build filename
		fnOut=pyfrp_misc_module.slashToFn(fnOut)+self.embryo.name+"_sim_"
		
		#Output
		print("Saving ", len(tvec), " images of simulation to ", fnOut)
		
		#Loop all tvec and build images
		j=0
		for i,t in enumerate(self.tvecSim):
			
			if t >= tvec[j]:
				j=j+1
				
				# Interpolate
				X,Y,img=self.computeInterpolatedSolutionToImg(self.vals[i],roi=roi,method=method,res=res)
				
				# If ROI is selected, fill everything outside of ROI with fillval
				if roi!=None:
					imgNew=fillVal*np.ones(img.shape)
					imgNew[roi.imgIdxX,roi.imgIdxY]=img[roi.imgIdxX,roi.imgIdxY]
					img=imgNew
					
				# Save img
				enum=(len(str(len(tvec)))-len(str(j)))*"0"+str(j)
				pyfrp_img_module.saveImg(img,fnOut+"t"+enum+".tif",scale=scale,enc=enc)
				
				if showProgress:
					print("Saved for t=", t ," to ", fnOut+"t"+enum+".tif")
			
		return True
		
		
	def visualize(self,cut=False,app=None):
		
		"""Visualizes simulation using VTK.
		
		Uses :py:class:`pyfrp.gui.pyfrp_gui_vtk.vtkSimVisualizer` to visualize the simulation. If ``cut=True``, will also
		add a cutter that allows to cut through the 3D simulation and display single slices, see also
		:py:class:`pyfrp.gui.pyfrp_gui_vtk.vtkSimVisualizerCutter`. Cutter only works for 3D simulations.
		
		If no application is specified, will launch a new ``QtWidgets.QApplication``. 
		
		.. note:: Simulation results must be saved in ``simulation`` object via ``saveSim=True``.
		
		Example:
		
		>>> emb.simulation.visualize(cut=False)
		
		.. image:: ../imgs/pyfrp_gui_vtk/vtkSimVisualizer.png
		
		and 
	
		>>> emb.simulation.visualize(cut=True)
		
		.. image:: ../imgs/pyfrp_gui_vtk/vtkSimVisualizerCutter.png
	
		Keyword Args:
			cut (bool): Visualize with cutter.
			app (QtWidgets.QApplication): Some application.
		
		"""
		
		# Check if running makes sense
		if len(self.vals)==0:
			printWarning("Embryo does not have saved simulation. Rerun simulation with saveSim=True.")
			return 
		
		if self.embryo.geometry.dim==2 and cut:
			printWarning("vtkSimVisualizerCutter does only work and make sense for 3D geometries. Falling back to vtkSimVisualizer.")
			cut=False
	
		# Import right module from gui submodule
		from pyfrp.gui import pyfrp_gui_vtk
		
		# Create QApplication
		if app==None:
			app = QtWidgets.QApplication([])
		
		if cut:
			widget=pyfrp_gui_vtk.vtkSimVisualizerCutter(self.embryo)
		else:
			widget=pyfrp_gui_vtk.vtkSimVisualizer(self.embryo)
		
		app.exec_()
	