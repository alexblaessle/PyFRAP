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

"""Simulaton module for PyFRAP toolbox. Handles simulating FRAP experiments and all necessary functions to do so, such as
	
	* Handling initial conditions.
	* Mimicing bleaching effects.
	* Experiment simulation.

"""

#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

#PDE Toolbox
from fipy import *

#Numpy/Scipy
import numpy as np
import scipy.interpolate as interp 
import scipy.ndimage.interpolation as ndi

#matplotlib
import matplotlib.pyplot as plt

#Misc
import time
import sys

#PyFRAP Modules
from . import pyfrp_plot_module 
from . import pyfrp_integration_module
from . import pyfrp_misc_module
from .pyfrp_term_module import *
from . import pyfrp_idx_module

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

def simulateReactDiff(simulation,signal=None,embCount=None,showProgress=True,debug=False):
	
	r"""Simulates reaction diffusion equation goverining FRAP experiment.
	
	Performs the following steps:
		
		* Resets ``simVecs`` of all ROIs.
		* If not generated yet, generates mesh (should never be the case!)
		* Initializes PDE with Neumann boundary conditions resulting in the problem:
		
		  .. math::
		     \partial_t c = D \nabla^2 c - k_1 c + k_2,
	
		  where :math:`k_1` is the degradation rate and :math:`k_2` the production rate.
		
		* Applies initial conditions defined in ``simulation.ICmode``.
		* Simulates FRAP experimment.
	
	Args: 
		simulation (pyfrp.subclasses.pyfrp_simulation.simulation): Simulation object.
	
	Keyword Args:
		signal (PyQt4.QtCore.pyqtSignal): PyQT signal to send progress to GUI.
		embCount (int): Counter of counter process if multiple datasets are analyzed. 
		debug (bool): Print final debugging messages and show debugging plots.
		showProgress (bool): Show simulation progress. 
		
	Returns: 
		pyfrp.subclasses.pyfrp_simulation.simulation: Updated simulation object.
	"""
	

	#Stepping and timescale
	timeStepDuration = simulation.tvecSim[1]-simulation.tvecSim[0]
	
	#Reset simulation vecs
	for r in simulation.embryo.ROIs:
		r.resetSimVec()
	
	#Empty list to put simulation values in
	vals=[]
	
	print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	print("Starting simulation")
	print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	
	startTimeTotal=time.clock()
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Mesh Generation
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	startTimeMesh=time.clock()
	
	if simulation.mesh.mesh==None:
		printWarning('No mesh has been generated yet!')
		a=input('Do you want to generate a mesh now?[Y/N]')
		if a=='Y':
			simulation.mesh.genMesh()
			print("Mesh created in", time.clock()-startTimeMesh)
		else:
			print('Cannot run simulation without mesh, will abort.')
			return simulation
		
		timeMesh=time.clock()-startTimeMesh
		
	print("Mesh created after", time.clock()-startTimeTotal)

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Initialization of PDE
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	#Create solution variable
	
	phi = CellVariable(name = "solution variable",mesh = simulation.mesh.mesh,value = 0.) 
	

	#Apply initial conditions
	if simulation.ICmode==0:
		phi = applyROIBasedICs(phi,simulation)
		
	elif simulation.ICmode==1:
		phi = applyRadialICs(phi,simulation,debug=debug)
	
	elif simulation.ICmode==2:
		phi=applyImperfectICs(phi,simulation,simulation.embryo.geometry.getCenter(),100.,simulation.embryo.sliceHeightPx)
		
	elif simulation.ICmode==3:
		phi=applyInterpolatedICs(phi,simulation,debug=False)
		
	elif simulation.ICmode==4:
		phi=applyIdealICs(phi,simulation,bleachedROI=simulation.bleachedROI,valOut=simulation.valOut)
		
	#Remember ICs
	simulation.IC=np.asarray(phi.value).copy()
	
	#Defining Type of equation
	eq = TransientTerm() == DiffusionTerm(coeff=simulation.D)+simulation.prod-simulation.degr*phi

	#Defining BCs
	#Note: BCs are Neumann boundaries by default 
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Calculating initial concentrations 
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	for r in simulation.embryo.ROIs:
		r.getSimConc(phi,append=True)
	
	if simulation.saveSim:
		vals.append(np.asarray(phi.value).copy())
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Solving PDE
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	#Keeping track of time
	startTimeSim=time.clock()
	
	avgTime=0
	stepTime=0
	
	#Choose solver
	if simulation.solver=="LU":
		mySolver = LinearLUSolver(iterations=simulation.iterations, tolerance=simulation.tolerance)
	elif simulation.solver=="PCG":
		mySolver = LinearPCGSolver(tolerance=simulation.tolerance,iterations=simulation.iterations)

	for step in range(simulation.stepsSim-1):
		
		#Compute timestep duration 
		timeStepDuration=simulation.tvecSim[step+1]-simulation.tvecSim[step]
		
		#Solve PDE in this Step
		stepStart=time.clock()
		eq.solve(var=phi,dt=timeStepDuration,solver=mySolver)
		stepTime=stepTime+(time.clock()-stepStart)
				
		#Compute concentration
		avgStart=time.clock()
		
		for r in simulation.embryo.ROIs:
			r.getSimConc(phi,append=True)
		
		avgTime=avgTime+(time.clock()-avgStart)
		
		#Save simulation array if necessary
		if simulation.saveSim:
			vals.append(np.asarray(phi.value).copy())
		
		#Print Progress
		if showProgress:
			currPerc=int(100*step/float(simulation.stepsSim))
			
			if signal==None:
				sys.stdout.write("\r%d%%" %currPerc)  
				sys.stdout.flush()
			else:	
				if embCount==None:
					signal.emit(currPerc)
				else:
					signal.emit(currPerc,embCount)
			

	print("Step time: ", stepTime, " in %:", stepTime/(time.clock()-startTimeSim)*100)
	print("Avg time: ", avgTime, " in %:", avgTime/(time.clock()-startTimeSim)*100)
	print("Simulation done after", time.clock()-startTimeTotal)
	
	#Save to simulation object only
	if simulation.saveSim:
		simulation.vals=list(vals)
	
	return simulation

def rerunReactDiff(simulation,signal=None,embCount=None,showProgress=True,debug=False):
	
	"""Reruns simulation by extracting values from ``simulation.vals``.
	
	Performs the following steps:
		
		* Resets ``simVecs`` of all ROIs.
		* Extracts values per ROI from ``simulation.vals``.
	
	.. note:: Only works if simulation has been run before with ``saveSim`` enabled.
	
	Args: 
		simulation (pyfrp.subclasses.pyfrp_simulation.simulation): Simulation object.
	
	Keyword Args:
		signal (PyQt4.QtCore.pyqtSignal): PyQT signal to send progress to GUI.
		embCount (int): Counter of counter process if multiple datasets are analyzed. 
		debug (bool): Print final debugging messages and show debugging plots.
		showProgress (bool): Show simulation progress. 
		
	Returns: 
		pyfrp.subclasses.pyfrp_simulation.simulation: Updated simulation object.
	"""
	
	#Check if can be rerun
	if len(simulation.vals)==0:
		printWarning("Values have not been saved for this simulation. Turn on saveSim to do that. Won't do anything for now.")
		return simulation
	
	#Reset simulation vecs
	for r in simulation.embryo.ROIs:
		r.resetSimVec()
	
	#Loop through vals
	for i,val in enumerate(simulation.vals):
		for r in simulation.embryo.ROIs:
			r.getSimConc(val,append=True)
		
		#Print Progress
		if showProgress:
			currPerc=int(100*i/float(len(simulation.vals)))
			
			if signal==None:
				sys.stdout.write("\r%d%%" %currPerc)  
				sys.stdout.flush()
			else:	
				if embCount==None:
					signal.emit(currPerc)
				else:
					signal.emit(currPerc,embCount)
		
	return simulation	
	

def applyROIBasedICs(phi,simulation):
	
	"""Applies ROI-based initial conditions.
	
	First sets concentration on all mesh nodes equal to `simulation.embryo.analysis.concRim`.
	Afterwards, mesh nodes get assigned the value of the first entry ``dataVec`` of 
	the ROI covering them. Note: If a mesh node is covered by two ROIs, will assign the value
	of the ROI that is last in embryo's ``ROIs`` list. See also 
	:py:func:`pyfrp.subclasses.pyfrp_simulation.setICMode`.
	
	Args:
		phi (fipy.CellVariable): PDE solution variable.
		simulation (pyfrp.subclasses.pyfrp_simulation.simulation): Simulation object.
	
	Returns:
		fipy.CellVariable: Updated solution variable.	
	
	"""
	
	phi.setValue(simulation.embryo.analysis.concRim)
	
	for r in simulation.embryo.ROIs:
		phi.value[r.meshIdx]=r.dataVec[0]
		
	
def applyIdealICs(phi,simulation,bleachedROI=None,valOut=None):
	
	"""Applies ideal initial conditions.
	
	That is, everything falling inside the bleached ROI in 
	x-y-direction will be set its initial dataVec value,
	everything else will be set equal to valOut.
	
	.. note:: The ``bleachedROI`` and valOut are often stored inside the simulation
	   object. If those two cannot be found, will try to find a ROI called *Bleached Square*
	   for the bleached ROI and set valOut to ``concRim``. If this again fails, will return
	   error.
	
	Args:
		phi (fipy.CellVariable): PDE solution variable.
		simulation (pyfrp.subclasses.pyfrp_simulation.simulation): Simulation object.
		
	Keyword Args:
		bleachedROI (pyfrp.subclasses.pyfrp_ROI.ROI): Bleached ROI.
		valOut (float): Value to be assigned outside of bleached ROI.
		
	Returns:
		fipy.CellVariable: Updated solution variable.	
	
	"""
	
	if bleachedROI==None:
		bleachedROI=simulation.embryo.getROIByName("Bleached Square")
		if bleachedROI==None:
			printError("No bleachedROI can be found in applyIdealICs.")
			return phi
		
	if valOut==None:
		if simulation.embryo.analysis.concRim!=None:
			valOut=simulation.embryo.analysis.concRim
		else:
			printError("No bleachedROI can be found in applyIdealICs.")
			return phi
	
	#Set all values in mesh to valOut
	phi.setValue(valOut)
	
	#Find indices of all nodes that lie inside bleachedROI in x-y-direction
	x,y,z=simulation.mesh.getCellCenters()
	ind=bleachedROI.checkXYInside(np.asarray(x),np.asarray(y))
	
	#Set these nodes equal to first dataVec entry
	phi.value[ind]=bleachedROI.dataVec[0]
	
	return phi
			
def applyRadialICs(phi,simulation,radSteps=15,debug=False):
	
	"""Applies radially averaged image data to solution variable as IC.
	
	.. note:: Will use ``embryo.geometry.center`` as center circle and the maximum 
	   distant pixel from the center as maximum radius.
	
	Args: 
		phi (fipy.CellVariable): PDE solution variable.
		simulation (pyfrp.subclasses.pyfrp_simulation.simulation): Simulation object.
	
	
		radSteps (int): Number of radial levels.
		debug (bool): Print debugging messages.
		
	Returns:
		fipy.CellVariable: Updated solution variable.
	
	"""
	
	#Adjust center so histogram works for 'quad'
	if 'quad' in list(simulation.embryo.analysis.process.keys()):
		center=[0,0]
	else:
		center=simulation.ICimg,simulation.embryo.geometry.getCenter()
	
	#Compute radial histogram of IC image
	maxR=pyfrp_img_module.dist(center,[simulation.ICimg.shape[0],simulation.ICimg.shape[0]])
	bins,binsMid,histY,binY=pyfrp_img_module.radialImgHist(simulation.ICimg,nbins=radSteps,byMean=True,maxR=maxR)
	
	#Set center to actual center of geometry
	center=simulation.ICimg,simulation.embryo.geometry.getCenter()
	
	#Apply value of most outer bin to all nodes
	phi.setValue(binY[-1])
	
	#Loop through all bins and apply values from outside to inside
	binsY=binsY.reverse()
	bins=bins.reverse()
	for i in range(len(binY)):
		
		phi.setValue(binY[i], where=(x-center[0])**2+(y-center[1])**2 < bins[i]**2)
			
		if debug:
			print("Applied concentration", binY[i], " to all nodes with radius <", bins[i]) 
	
	return phi
		

def applyInterpolatedICs(phi,simulation,matchWithMaster=True,debug=False,fixNeg=True,fillICWithConcRim=True):
	
	"""Interpolates initial conditions onto mesh.
	
	Uses a bivarariate spline interpolation (http://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.interpolate.RectBivariateSpline.html)
	to generate an interpolation function of the IC image. Then applies interpolated values to solution variable ``phi`` if mesh nodes are inside 
	image and masterROI. If not, will apply rim concentration.
	
	.. note:: If no rim concentration has been calculated (for example through running the data analysis)
	   applyInterpolatedICs will try to compute ``concRim`` by itself. For this it will take the mean concentration outside of bleached square but inside
	   ``masterROI``. 
	   
	.. note:: The bleached square used here is not defined as a ``ROI`` object here, but rather through the properties 
	   ``embryo.sideLengthBleachedPx`` and ``embryo.offsetBleachedPx``. This might change in future versions.
	
	Args: 
		phi (fipy.CellVariable): PDE solution variable.
		simulation (pyfrp.subclasses.pyfrp_simulation.simulation): Simulation object.
	
	Keyword Args:
		matchWithMaster (bool): Match interpolation indices with ``masterROI`` indices.
		debug (bool): Print debugging messages.
	
	Returns:
		fipy.CellVariable: Updated solution variable.
	"""
	
	#Get image resolution and center of geometry
	res=simulation.ICimg.shape[0]
	center=simulation.embryo.geometry.getCenter()
	
	#Define x/y coordinates of interpolation
	if 'quad' in list(simulation.embryo.analysis.process.keys()):
		#Shift everything by center to fit with the mesh
		xInt = np.arange(center[0]+1, center[0]+res+1, 1)
		yInt = np.arange(center[1]+1, center[1]+res+1, 1)		
	else:
		xInt = np.arange(1, res+1, 1)
		yInt = np.arange(1, res+1, 1)
	
	#Getting cell centers
	x,y,z=simulation.mesh.getCellCenters()
	
	#Finding outer rim concentration
	if simulation.embryo.analysis.concRim==None:
		printWarning('concRim is not analyzed yet. Will use concentration outside of bleached region as approximation')
		
		#Grab offset and sidelength
		offset=simulation.embryo.offsetBleachedPx
		sidelength=simulation.embryo.sideLengthBleachedPx
		
		#Get indices outside of bleached square but inside masterROI
		indXSqu,indYSqu=pyfrp_idx_module.getSquareIdxImg(simulation.embryo.offsetBleachedPx,simulation.embryo.sideLengthBleachedPx,simulation.embryo.dataResPx)
		
		masterROI=simulation.embryo.getMasterROI()
		
		indX=pyfrp_misc_module.complValsSimple(masterROI.imgIdxX,indXSqu)
		indY=pyfrp_misc_module.complValsSimple(masterROI.imgIdxX,indYSqu)
		
		if 'quad' in list(simulation.embryo.analysis.process.keys()):
			img=pyfrp_img_module.unflipQuad(np.flipud(simulation.ICimg))
		else:
			img=simulation.ICimg
		
		concRim=pyfrp_img_module.meanConc(img[indX,indY])
		
		print('Approximate concRim = ', concRim)
		
	else:	
		concRim=simulation.embryo.analysis.concRim
	
	if fillICWithConcRim:
		masterROI=simulation.embryo.getMasterROI()
		ICimg=concRim*np.ones(simulation.ICimg.shape)
		ICimg[masterROI.imgIdxX,masterROI.imgIdxY]=simulation.ICimg[masterROI.imgIdxX,masterROI.imgIdxY]	
	else:	
		ICimg=simulation.ICimg.copy()
	
	#Generate interpolation function
	f=interp.RectBivariateSpline(xInt, yInt, ICimg.T, bbox=[None, None, None, None], kx=3, ky=3, s=0)
	
	#Set all values of solution variable to concRim
	phi.setValue(concRim)
	
	#Get Offset of image and check which nodes are inside image
	if 'quad' in list(simulation.embryo.analysis.process.keys()):
		offset=[simulation.embryo.dataResPx/2,simulation.embryo.dataResPx/2]
		ins=pyfrp_idx_module.checkInsideImg(x,y,simulation.embryo.dataResPx/2,offset=offset)
	else:
		offset=[0,0]
		ins=pyfrp_idx_module.checkInsideImg(x,y,simulation.embryo.dataResPx,offset=offset)
		
	#Convert into indices
	ind=np.arange(len(x))
	ind=ind[np.where(ins)[0]]

	"""NOTE:  I think we need to match here indices inside image with the one of master ROI, so we don't apply
	values outside of masterROI (generally just background) to nodes that lie INSIDE image, but OUTSIDE of masterROI.
	"""
	
	if matchWithMaster:
		
		masterROI=simulation.embryo.getMasterROI()
	
		xnew=np.asarray(x)[ind]
		ynew=np.asarray(y)[ind]
				
		ins=masterROI.checkXYInside(xnew,ynew)
		
		ind=np.asarray(ind)
		ind=ind[np.where(ins)[0]]
		ind=list(ind)
		
	#Apply interpolation
	try:
		phi.value[ind]=f.ev(x[ind],y[ind])
	except IndexError:
		if debug:
			printNote("Changed index array to nparray b/c of IndexError.")
		ind=np.array(ind)
		phi.value[ind]=f.ev(x[ind],y[ind])
	
	#Fix negative values if selected
	if fixNeg:
		phi=fixNegValues(phi)
	
	return phi
	
def fixNegValues(phi,minVal=None):
	
	"""Fixes negative values in solution variable. 
	
	Interpolation sometimes returns negative values if gradients are really steep. Will
	apply ``minVal`` to such nodes. 
	
	If ``minVal==None``, will take smallest non-negative value of solution value.
	
	
	
	"""
	
	if minVal==None:
		minVal=min(phi.value[np.where(phi.value>=0)[0]])
	
	phi.value[np.where(phi.value<0)[0]]=minVal
	
	return phi
	
def sigmoidBleachingFct(x,y,z,center,rJump,sliceHeight,maxVal=1.,maxMinValPerc=0.25,minMinValPerc=0.25,rate=0.1):
	
	r"""Generates sigmoid scaling function for imperfect bleaching at
	coordinates x/y/z.
	
	The idea behind the sigmoid function is:
	
		* Through scattering and other effects, the bleached window becomes blurry in larger depths, resulting
		  in a radial sigmoid scaling function around ``center``. 
		* Similarly, bleaching intensity increases with depth. Thus, a linear term controls the values close 
		  to ``center`` of the sigmoid function. Bleaching closer to the laser than the imaged height will 
		  be rendered stronger, while bleaching effects below will be decreased by *bumping up* the 
		  bleached window. However, only until some threshold is reached.
		  
	The sigmoid function is given by:	  
	
	.. math:: s(r,z) = v_{\mathrm{min}}(z)+(v_{\mathrm{max}}-v_{\mathrm{min}}(z))\frac{1}{1+\exp(-\rho(r-r_\mathrm{Jump}))},
	
	where :math:`\rho` is the sigmoid slope given by ``rate``, :math:`r_\mathrm{Jump}` is the radius from ``center``
	at which sigmoid function has its **jump**, given by ``rJump`` and :math:`r` the radius of coordinate ``[x,y]`` from
	``center``.
	
	:math:`v_{\mathrm{min}}(z)` is a linear function describing how strong the bleaching is dependent on the 
	depth :math:`z` given by
	
	.. math:: v_{\mathrm{min}}(z) = \frac{v_{\mathrm{max}} - v_{\mathrm{max. bleach}}}{h_s} z +  v_{\mathrm{max. bleach}},
	
	where :math:`v_{\mathrm{max}}` is the value of the sigmoid function far from ``center``, :math:`v_{\mathrm{max. bleach}}`
	is the strongest expected bleaching value (mostly at :math:`z=0`) and :math:`h_s` is the height of the imaging slice, given 
	by ``sliceHeight``. 
	The maximum rate of bleaching :math:`v_{\mathrm{max. bleach}}` is computed by:
	
	.. math:: v_{\mathrm{max. bleach}} = (1-p_{\mathrm{max. bleach}})v_{\mathrm{max}},
	
	where :math:`p_{\mathrm{max. bleach}}` is the percentage of maximum expected bleaching compared to the values in the imaging 
	height, given by ``maxMinValPerc``. That is, by how much has the laser power already decreased on its way from entry point of 
	the sample to the imaging height.
	
	For sample depths deeper than the imaging height, bleaching is expected to be decrease in intensity, thus the bleached area
	is getting **bumped up**. To avoid bumping the bleached area by too much, eventually even resulting in the bleached
	area having a higher concentration than the area outside, the sigmoid function has a cut-off: If values of :math:`s(r,z)` pass 
	
	.. math:: v_{\mathrm{min. bleach}} = (1+p_{\mathrm{min. bleach}})v_{\mathrm{max}},
	
	where :math:`p_{\mathrm{min. bleach}}` is the percentage of bleaching to cut-off, then we set
	
	.. math:: s(r,z) = v_{\mathrm{min. bleach}},
	
	ultimately resulting in a scaling function given by
	
	.. math:: s(r,z) = \left\{\begin{array}{cc} 
	   v_{\mathrm{min}}(z)+(v_{\mathrm{max}}-v_{\mathrm{min}}(z))\frac{1}{1+\exp(-\rho(r-r_\mathrm{Jump}))} & \mbox{ if } s(r,z) <= v_{\mathrm{min. bleach}} , \\
	   v_{\mathrm{min. bleach}} & \mbox{ else }
	   \end{array}
	   \right.
	   
	.. image:: ../imgs/pyfrp_sim_module/sigmoidFct.png
	
	Args:
		x (numpy.ndarray): x-coordinates.
		y (numpy.ndarray): y-coordinates.
		z (numpy.ndarray): z-coordinates.
		center (list): Center of bleaching. 
		rJump (float): Radius from center where sigmoid jump is expected.
		sliceHeight (float): Height at which dataset was recorded. 
		
	Keyword Args:	
		maxVal (float): Value of sigmoid function outside of bleached region.
		maxMinValPerc (float): Percentage of maximum bleaching intensity.
		minMinValPerc (float): Percentage of minimum bleaching intensity. 
		rate (float): Rate at which sigmoid increases.
			
	Returns:
		numpy.ndarray: 
		
	"""
	
	#Calculate linear equation describing how strong bleaching effect decreases as z increases
	maxMinVal=(1-maxMinValPerc)*maxVal
	m=(maxVal-maxMinVal)/sliceHeight
	b=maxMinVal
	minVal=m*z+b
	
	#Compute distance from center for each point
	r=np.sqrt((x-center[0])**2+(y-center[1])**2)
	
	#Compute sigmoid function
	sigm=minVal+(maxVal-minVal)/(1+np.exp(-rate*(r-rJump)))
	
	#Compute cut-off, so we do not amplifiy the bleached region
	minMinVal=(1+minMinValPerc)*maxVal
	sigm[np.where(sigm>minMinVal)]=minMinVal
	
	return sigm,r
	
def applyImperfectICs(phi,simulation,center,rJump,sliceHeight,maxVal=1.,maxMinValPerc=0.25,minMinValPerc=None,rate=0.1,matchWithMaster=True,debug=False):
	
	"""Mimic imperfect bleaching through cone approximation, return phi.
	
	 .. warning:: Not working in current version. Will be integrated in further versions again.
	
	"""
	
	phi = applyInterpolatedICs(phi,simulation,matchWithMaster=matchWithMaster,debug=debug)
	
	#Get cell coordinates
	x,y,z=simulation.mesh.getCellCenters()
	x=np.asarray(x)
	y=np.asarray(y)
	z=np.asarray(z)
	
	#Compute by how much 
	if minMinValPerc==None:
		r=np.sqrt((x-center[0])**2+(y-center[1])**2)
		inVal=np.mean(phi.value[np.where(r<rJump)])
		outVal=np.mean(phi.value[np.where(r>=rJump)])
		
		minMinValPerc=inVal/outVal
		

		
	#Compute sigmoid function
	sigm,r = sigmoidBleachingFct(x,y,z,center,rJump,sliceHeight,maxVal=1.,maxMinValPerc=maxMinValPerc,minMinValPerc=minMinValPerc,rate=rate)
	
	#Multiplicate solution variable with sigmoid function
	phi.value = phi.value * sigm

	return phi
	
	



		
	

