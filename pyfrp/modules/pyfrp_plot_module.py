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

"""Plotting module for PyFRAP toolbox. 

Contains functions and classes that are often used by PyFRAP toolbox and simplify plot creation and management.

"""

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#numpy
import numpy as np
import scipy.interpolate

#Plotting
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.patches as ptc

#Misc
import sys
import os

#PyFRAP
import pyfrp_img_module
from pyfrp_term_module import *
from pyfrp.modules import pyfrp_idx_module

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

class FRAPBoundarySelector():
	
	"""Simple GUI widget to select circular FRAP boundaries.
	
	Has useful center marker that is activatable that helps finding the 
	center of the image.
	
	Mouse Input:
		
		* Left: Set center.
		* Right: Set Radius.
		* Middle: Activate center marker.
		
	Keyboard Input:
	
		* Left Arrow: Move center to the left.
		* Right Arrow: Move center to the right.
		* Up Arrow: Move center upwards.
		* Down Arrow: Move center downwards.
		* Control + Up Arrow: Increase circle radius.
		* Control + Down Arrow: Decrease circle radius.
	
	.. note:: If ``embryo`` is given at initiation, will use first image specified
	   in embryo's ``fileList`` as background image.
	   
	Example Usage:
	
	>>> sel=FRAPBoundarySelector(fn="path/to/img/file")

	Use mouse/keyboard to define circular boundary.
	
	>>> center,radius=sel.getResults()
	
	Keyword Args:
		embryo (pyfrp.subclasses.pyfrp_embryo.embryo): PyFRAP embryo instance.
		fn (str): Filepath to image file taken for boundary selection.
		
	"""
	
	def __init__(self,embryo=None,fn=None,img=None):
		
		#Passing embryo to class
		self.embryo=embryo	
		self.fn=fn
		
		#Img
		self.img=img
		
		#Plot resolution
		self.dpi = 100
	
		#Some bookkeeping variables
		self.radiusPt=[]
		self.centerPt=[]
		self.centerMarker=[]
		
		#Results
		self.radius=None
		self.center=None
		
		#Creating figure and canvas
		self.createCanvas()
		
		#Close everthing down if faulty input
		self.checkInput()
			
		#Plot image if existent
		self.showFirstDataImg()
		plt.show()
	
	def checkInput(self):
		
		"""Checks if at least one of the two keyword arguments, ``embryo`` or ``fn``, is given.
		
		If not, prints error message and closes down widget.
		"""
		
		if self.fn==None and self.embryo==None and self.img==None:
			printError("No Embryo, image or fn defined. Going to exit.")
			plt.close(self.fig)
			return
	
	def createCanvas(self):
		
		"""Creates figure and canvas used for plotting."""
			
		h=500/self.dpi
		v=500/self.dpi
		
		self.fig = plt.figure()
		self.fig.show()
		self.fig.set_size_inches(h,v,forward=True)
		
		self.canvas=self.fig.canvas
		
		self.ax = self.fig.add_subplot(111)
		
		if self.embryo!=None:
			self.ax.set_xlim([0, self.embryo.dataResPx])
			self.ax.set_ylim([0, self.embryo.dataResPx])
			
		#Connect to mouseclick event
		self.fig.canvas.mpl_connect('close_event', self.closeFigure)
		self.canvas.mpl_connect('button_press_event', self.getMouseInput)
		self.canvas.mpl_connect('key_press_event', self.keyPressed)
		
		self.canvas.draw()
		
		return 
	
	def keyPressed(self,event):
		
		"""Directs all key press events to the respective functions."""
		
		if event.key=='left':
			self.moveLeft()
		elif event.key=='right':
			self.moveRight()
		elif event.key=='up':
			self.moveUp()
		elif event.key=='down':
			self.moveDown()
		elif event.key=='ctrl+up':
			self.increaseRadius()
		elif event.key=='ctrl+down':
			self.decreaseRadius()
			
	def moveLeft(self):
		
		"""Moves center 1 px to the left."""
		
		if self.center!=None:
			self.center=[self.center[0]-1,self.center[1]]
			self.drawCenter()
		
		
	def moveRight(self):
		
		"""Moves center 1 px to the right."""
		
		if self.center!=None:
			self.center=[self.center[0]+1,self.center[1]]
			self.redraw()
	
	def moveUp(self):
		
		"""Moves center 1 px up."""
		
		if self.center!=None:
			self.center=[self.center[0],self.center[1]+1]
			self.redraw()
	
	def moveDown(self):
		
		"""Moves center 1 px down."""
		
		if self.center!=None:
			self.center=[self.center[0],self.center[1]-1]
			self.redraw()
	
	def redraw(self):
		
		"""Redraws both center and radius if available."""
		
		if self.center!=None:
			self.drawCenter()
		if self.radius!=None:
			self.drawRadius()
			
	def increaseRadius(self):
		
		"""Increases radius by 1 px."""
		
		if self.radius!=None:
			self.radius=self.radius+1
			self.redraw()
			
	def decreaseRadius(self):
		
		"""Decreases radius by 1 px."""
		
		if self.radius!=None:
			self.radius=self.radius-1
			self.redraw()
			
	def closeFigure(self,event):
		
		"""Returns center and radius at close ``event``."""
		
		return self.center,self.radius
		
	def getEmbryo(self):	
		
		"""Returns ``embryo`` object if given.``"""
		
		return self.embryo
	
	def getResults(self):
		
		"""Returns center and radius."""
		
		return self.center,self.radius
	
	def drawCenterMarker(self):
		
		"""Draws a yellow marker in center of the image, making it 
		easier to find image center when selecting center of boundary."""
		
		centerImg=[self.img.shape[0]/2.,self.img.shape[1]/2.]
		
		if len(self.centerMarker)>0:
			self.clearCenterMarker()
		else:	
			pt=ptc.Circle(centerImg,radius=3,fill=True,color='y')
			self.centerMarker.append(self.ax.add_patch(pt))
			
		self.fig.canvas.draw()
		
		return self.centerMarker
		
	def clearCenterMarker(self):
		
		"""Removes center maker from canvas."""
	
		for pt in self.centerMarker:
			pt.remove()
			
		self.fig.canvas.draw()	
			
		self.centerMarker=[]
	
	def drawCenter(self):
		
		"""Draws a red marker at selected center on canvas."""
		
		if len(self.centerPt)>0:
			self.clearCenter()
	
		pt=ptc.Circle(self.center,radius=3,fill=True,color='r')
		self.centerPt.append(self.ax.add_patch(pt))
		
		self.fig.canvas.draw()
		
		return self.centerPt
	
	def clearCenter(self):
	
		"""Removes center marker from canvas."""
		
		for pt in self.centerPt:
			pt.remove()
			
		self.centerPt=[]
		
		self.fig.canvas.draw()
		
	def drawRadius(self):
		
		"""Draws a red circle around selected center with selected radius on canvas."""
		
		if len(self.radiusPt)>0:
			self.clearRadius()
			
		pt=ptc.Circle(self.center,radius=self.radius,fill=False,color='r')
		self.radiusPt.append(self.ax.add_patch(pt))
		
		self.fig.canvas.draw()
		
		return self.radiusPt
	
	def clearRadius(self):
		
		"""Removes circle from canvas."""
		
		for pt in self.radiusPt:
			pt.remove()
			
		self.radiusPt=[]
		
		self.fig.canvas.draw()
		
		return self.radiusPt
		
		
	def showFirstDataImg(self):
		
		"""Shows either first data image defined in ``embryo.fileList`` or 
		image specified by ``fn``.
		
		.. note:: If both are given, will use the embryo option.
		
		"""
		
		if self.embryo!=None:
			self.embryo.updateFileList()
		
			fnImg=self.embryo.fnDatafolder
			if fnImg[-1]!='/':
				fnImg=fnImg+'/'
			fnImg=fnImg+self.embryo.fileList[0]
			
			self.img=pyfrp_img_module.loadImg(fnImg,self.embryo.dataEnc)
		
		elif self.fn!=None:
			self.img=pyfrp_img_module.loadImg(self.fn,'uint16')
			self.ax.set_xlim([0, self.img.shape[0]])
			self.ax.set_ylim([0, self.img.shape[1]])
		
		
		
		self.showImg(self.img)
	
	def showImg(self,img):
		
		"""Shows image on canvas.
		
		Args:
			img (numpy.ndarray): Image to be shown.
		
		"""
		
		self.ax.imshow(img)
		self.fig.canvas.draw()
		
		return self.canvas
	
	def computeRadiusFromCoordinate(self,x,y):
		
		"""Computes radius from given cordinate ``(x,y)``.
		"""
		
		return np.sqrt((x-self.center[0])**2+(y-self.center[1])**2)
		
	def getMouseInput(self,event):
		
		"""Directs mouse input to the right actions."""
		
		#Check if click in axes
		if event.xdata==None:
			return
		
		#Left click to define center
		if event.button==1:
			
			self.center=[event.xdata,event.ydata]
			
			self.drawCenter()
			
			if len(self.radiusPt)>0:
				self.drawRadius()
				
		#Right click to define radius
		elif event.button==3:
			
			if len(self.centerPt)>0:
			
				self.radius=self.computeRadiusFromCoordinate(event.xdata,event.ydata)
				self.drawRadius()
			
		#Middle click to activate center marker	
		if event.button==2:
			
			self.drawCenterMarker()
				
		self.fig.canvas.draw()
		
		return
	
def makeSubplot(size,titles=None,tight=False,sup=None,proj=None,fig=None,show=True):
	
	"""Generates matplotlib figure with (x,y) subplots.
	
	.. note:: List of ``titles`` needs to be same size as desired number of axes.
	   Otherwise will turn off titles.
	
	.. note:: List of ``proj`` needs to be same size as desired number of axes.
	   Otherwise will turn off projections.
	   
	Example:
	
	>>> makeSubplot([2,2],titles=["Axes 1", "Axes 2", "Axes 3", "Axes 4"],proj=[None,None,'3d',None])
	
	will result in
	
	.. image:: ../imgs/pyfrp_plot_module/makeSubplot.png
	
	Args:
		size (list): Size of subplot arrangement.
	
	Keyword Args:
		titles (list): List of axes titles.
		tight (bool): Use tight layout.
		sup (str): Figure title.
		proj (list): List of projections.
		fig (matplotlib.figure): Figure used for axes.
		show (bool): Show figure right away.
	
	Returns:
		tuple: Tuple containing:
			
			* fig (matplotlib.figure): Figure.
			* axes (list): List of Matplotlib axes.
	
	"""
	
	#How many axes need to be created
	n_ax=size[0]*size[1]
	
	if proj==None:
		proj=n_ax*[None]
	
	#Creating figure
	if fig==None:
		fig=plt.figure()
		if show:
			fig.show()
	fig.set_tight_layout(tight)
	
	#Suptitle
	if sup!=None:
		fig.suptitle(sup)
	
	#Add axes
	axes=[]
	for i in range(n_ax):
		
		try:
			if proj[i]!=None:
				ax=fig.add_subplot(size[0],size[1],i+1,projection=proj[i])	
			else:
				ax=fig.add_subplot(size[0],size[1],i+1)
		except IndexError:
			printWarning("Axes " + str(i) + "does not have a projection specified. Will create normal 2D plot.")
			ax=fig.add_subplot(size[0],size[1],i+1)
				
		axes.append(ax)
		
		#Print titles 
		if titles!=None:
			try:
				ax.set_title(titles[i])
			except IndexError:
				printWarning("Axes " + str(i) + " does not have a title specified.")
	#Draw
	plt.draw()
	
	#Return axes handle
	return fig,axes

def makeGeometryPlot(titles=None,tight=False,sup=None,fig=None,show=True,unit="px"):
	
	"""Generates matplotlib figure  and single axes optimized for geometry plots.
	
	See also :py:func:`pyfrp.modules.pyfrp_plot_module.makeSubplot`.
	
	Keyword Args:
		titles (list): List of axes titles.
		tight (bool): Use tight layout.
		sup (str): Figure title.
		fig (matplotlib.figure): Figure used for axes.
		show (bool): Show figure right away.
		unit (str): Unit displayed in axes label.
	
	Returns:
		tuple: Tuple containing:
			
			* fig (matplotlib.figure): Figure.
			* axes (list): List of Matplotlib axes.
	
	"""
	
	fig,axes=makeSubplot([1,1],proj=['3d'],titles=titles,tight=tight,sup=sup,fig=fig,show=show)
	
	if unit==None:
		axes[0].set_xlabel("x")
		axes[0].set_ylabel("y")
		axes[0].set_zlabel("z")
	else:	
		axes[0].set_xlabel("x ("+unit+")")
		axes[0].set_ylabel("y ("+unit+")")
		axes[0].set_zlabel("z ("+unit+")")
	
	return fig,axes
	
	

def adjustImshowRange(axes,vmin=None,vmax=None):
	
	"""Adjust display range of ``matplotlib.pyplot.imshow`` plots in 
	list of axes.
	
	Finds first image artist in each axes in ``axes`` list and then
	sets display range to ``[vmin,vmax]``.
	
	Args:
		axes (list): List of matplotlib axes.
		
	Keyword Args:
		vmin (float): Minimum value of display range.
		vmax (float): Maximum value of display range.
		
	Returns:
		 list: Updated list of matplotlib axes.
	
	"""
	
	#Loop through axes
	for ax in axes:
		
		#Grab plot
		implot=findArtist(ax,"Image")
		
		#Rescale
		if implot!=None:
			implot.set_clim(vmin,vmax)
		
	#Draw
	redraw(ax)
		
	return axes

def findArtist(ax,key):
	
	"""Finds ``matplotlib.artist`` which name contains ``key``.
	
	.. note:: Will stop searching after first artist is found.
	
	Will return ``None`` if no artist can be found.
	
	Args:
		ax (matplotlib.axes): Matplotlib axes.
		key (str): Key used for search.
		
	Returns:
		matplotlib.artist: Matplotlib artist.
	
	"""
	
	c=ax.get_children()
	
	for x in c:
		if key in str(x):
			return x
			break
		
	return None	

def redraw(ax):
	
	"""Redraws axes's figure's canvas.
	
	Makes sure that current axes content is visible.
	
	Args:
		ax (matplotlib.axes): Matplotlib axes.
		
	Returns:
		matplotlib.axes: Matplotlib axes 
	
	"""
			
	ax.get_figure().canvas.draw()
	
	return ax

def plotTS(xvec,yvec,label='',title='',sup='',ax=None,color=None,linewidth=1,legend=True,linestyle='-',show=True,alpha=1.,legLoc=-1):
	
	"""Plot timeseries all-in-one function.
		
	Args:
		xvec (numpy.ndarray): x-data to be plotted.
		yvec (numpy.ndarray): y-data to be plotted.
	
	Keyword Args:
		ax (matplotlib.axes): Matplotlib axes used for plotting. If not specified, will generate new one.
		color (str): Color of plot.
		linestyle (str): Linestyle of plot.
		linewidth (float): Linewidth of plot.
		legend (bool): Show legend.
		sup (str): Figure title.
		title (str): Axes title.
		label (str): Label for legend.
		show (bool): Show figure.
		alpha (float): Transparency of line.
		legLoc (int): Location of legend.
	
	Returns:
		matplotlib.axes: Axes used for plotting.
	
	"""
	
	if len(xvec)!=len(yvec):
		printWarning('len(xvec) != len (yvec). This could be due to incomplete simulation/analysis/pinning. Will not plot.')
		return None
	
	if ax==None:
		fig,axes = makeSubplot([1,1],titles=[title],sup=sup,tight=False,show=show)
		ax=axes[0]
	else:
		ax.set_title(title)
		
	ax.plot(xvec,yvec,color=color,label=label,linestyle=linestyle,alpha=alpha)
	
	if legend:
		if legLoc==-1:
			ax.get_figure().subplots_adjust(right=0.7)
			ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		else:
			ax.legend(loc=legLoc)
		
	redraw(ax)
	
	return ax

def plotSolutionVariable(x,y,val,ax=None,vmin=None,vmax=None,nlevels=25,colorbar=True,plane='xy',zs=None,zdir=None,title="Solution Variable",sup="",dThresh=None,nPts=1000,mode='normal',typ='contour'):
		
	"""Plots simulation solution variable as 2D contour plot or 3D surface plot.
	
	.. note:: If no ``ax`` is given, will create new one.
	
	.. note:: ``x`` and ``y`` do not necessarily have to be coordinates in x/y-direction, but rather correspond to 
	   the two directions defined in ``plane``. 
	
	``plane`` variable controls in which plane the solution variable is supposed to be plotted. 
	Acceptable input variables are ``"xy","xz","yz"``. See also
	:py:func:`pyfrp.subclasses.pyfrp_ROI.ROI.getMaxExtendPlane`.
	
	``mode`` controls which contour/surface function is used:
	
		* ``normal``: Will create rectangular grid and use ``scipy.interpolate.gridddata`` and interpolate solution onto it,
		  then plot using ``matplotlib.pyplot.contourf``, see also http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.contourf and
		  http://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.interpolate.griddata.html#scipy.interpolate.griddata .
		
		* ``tri``: Will plot irregular grid using ``matplotlib.pyplot.tricontourf``. See also 
		  http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.tricontourf .
	
	``typ`` controls which plot type is used:
	
		* ``contour``: Produces contour plots using either ``matplotlib.pyplot.contourf`` or ``matplotlib.pyplot.tricontourf``.
		
		* ``surface``: Produces contour plots using either ``matplotlib.Axed3D.plot_surface`` or ``matplotlib.pyplot.plot_trisurf``.
	
	.. warning:: ``matplotlib.pyplot.tricontourf`` has problems when ``val`` only is in a single level of contour plot.
		To avoid this, we currently add some noise in this case just to make it plottable. This is not the most elegant
		solution. (only in case of 3D plotting)
	
	Args:
		x (numpy.ndarray): x-coordinates.
		y (numpy.ndarray): y-coordinates.
		val (numpy.ndarray): Solution variable values.
		
	Keyword Args:
		ax (matplotlib.axes): Axes used for plotting.
		vmin (float): Minimum value displayed in contour plot.
		vmax (float): Maximum value displayed in contour plot.
		nlevels (int): Number of contour levels.
		colorbar (bool): Display color bar.
		plane (str): Plane in which solution variable is supposed to be plotted.
		zs (float): In case of a 3D plot, height in direction zdir where to put contour.
		zdir (str): Orthogonal direction to plane.
		nPts (int): Number of points used for interpolating (only if ``mode=normal``).
		mode (str): Which contour function to use.
		title (str): Title of plot.
		typ (str): Type of plot.
		
	Returns:
		matplotlib.axes: Axes used for plotting.
	
	"""
	
	#Check of entered plot type makes sense
	if typ not in ['contour','surface']:
		printError("Unknown plot type "+ typ)
		return ax
	
	#Make axes if necessary
	if ax==None:
		if (zs!=None and zdir!=None) or typ=='surface':
			fig,axes = makeSubplot([1,1],titles=[title],sup=sup,proj=['3d'])
		else:	
			fig,axes = makeSubplot([1,1],titles=[title],sup=sup)
		
		ax=axes[0]
	else:
		ax.set_title(title)
		
	
	#vmin/vmax/levels
	vmin,vmax,levels=makeFittingLevels(vmin,vmax,val,nlevels)
	
	#Make grid
	grid = np.meshgrid(np.linspace(min(x), max(x), nPts),np.linspace(min(y), max(y), nPts))
	
	#Interpolate
	xy=np.vstack((x,y))
	valPlot = scipy.interpolate.griddata(xy.T, val, tuple(grid), 'linear')
	
	#If dThresh is given, filter all nodes that are further apart from their neighbor than dThresh
	if dThresh!=None:
		idxs=pyfrp_idx_module.maskMeshByDistance(x,y,dThresh,grid)
		valPlot[idxs]=np.nan
		
	#Stupid fix for tricontourf 3D problem
	if zs!=None and zdir!=None and mode=='tri':
		if min(val)==max(val):
			
			"""IDEA:
			Add a little bit of noise to make it work. Not best solution ever. Needs to be changed.
			"""

			add=(1.01*vmax-vmin)/(nlevels)*np.random.randn(np.shape(val)[0])
			val=val+add
	
	#Plot 
	if mode=='tri':
		if typ=='surface':
			solPlot=ax.plot_trisurf(x,y,val,cmap='jet',vmin=vmin,vmax=vmax)
		elif typ=='contour':
			solPlot=ax.tricontourf(x,y,val,vmin=vmin,vmax=vmax,levels=levels,offset=zs,zdir=zdir,extend='both')
	else:
		
		if typ=='surface':
			solPlot=ax.plot_surface(grid[0],grid[1],valPlot,cmap='jet',vmin=vmin,vmax=vmax)
		elif typ=='contour':
			solPlot=ax.contourf(grid[0],grid[1],valPlot,vmin=vmin,vmax=vmax,levels=levels,offset=zs,zdir=zdir)	
	
	#Label
	if len(plane)!=2:
		printError("Don't understand plane="+plane+". Will not plot.")
		return ax	
	else:
		ax.set_xlabel(plane[0])
		ax.set_ylabel(plane[1])

	ax.autoscale(enable=True, axis='both', tight=True)
	
	if colorbar:
		cb=plt.colorbar(solPlot,orientation='horizontal',pad=0.05)
	
	#plt.draw()
	ax.get_figure().canvas.draw()
	
	return ax	

def makeFittingLevels(vmin,vmax,val,nlevels,buff=0.01):
	
	"""Generates array with fitting levels for contour plots.
	
	.. note:: If ``vmin=None`` or ``vmax=None``, will pick the minimum/maximum 
	   value of ``val``.
	
	Args:
		val (numpy.ndarray): Array to be plotted.
		vmin (float): Minimum value displayed in contour plot.
		vmax (float): Maximum value displayed in contour plot.
		nlevels (int): Number of contour levels.
		
	Keyword Args:	
		buff (float): Percentage buffer to be added to both sides of the array.
		
	Returns:
		tuple: Tuple containing:
		
			* vmin (float): Minimum value displayed in contour plot.
			* vmax (float): Maximum value displayed in contour plot.
			* levels (numpy.ndarray): Level array to be handed to contour function.
	
	"""
	
	
	if vmin==None:
		vmin=min(val)
	if vmax==None:
		vmax=max(val)
		
	levels=np.linspace((1-buff)*vmin,(1+buff)*vmax,nlevels)
	
	return vmin,vmax,levels
	
def set3DAxesEqual(ax):

	"""Make axes of 3D plot have equal scale.
	
	This is one possible solution to Matplotlib's
	ax.set_aspect('equal') and ax.axis('equal') not working for 3D.
	
	Modified from http://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to .
	
	Args:
		ax (matplotlib.axes): A matplotlib axes.
		
	Returns:
		matplotlib.axes: Modified matplotlib axes.
		
	"""

	xLimits = ax.get_xlim3d()
	yLimits = ax.get_ylim3d()
	zLimits = ax.get_zlim3d()

	xRange = abs(xLimits[1] - xLimits[0])
	xMiddle = np.mean(xLimits)
	yRange = abs(yLimits[1] - yLimits[0])
	yMiddle = np.mean(yLimits)
	zRange = abs(zLimits[1] - zLimits[0])
	zMiddle = np.mean(zLimits)

	# The plot bounding box is a sphere in the sense of the infinity
	# norm, hence I call half the max range the plot radius.
	plotRadius = 0.5*max([xRange, yRange, zRange])

	ax.set_xlim3d([xMiddle - plotRadius, xMiddle + plotRadius])
	ax.set_ylim3d([yMiddle - plotRadius, yMiddle + plotRadius])
	ax.set_zlim3d([zMiddle - plotRadius, zMiddle + plotRadius])
	
	return ax

def getPubParms(fontSize=10,labelFontSize=10,tickFontSize=10,legendFontSize=10):
	
	"""Returns dictionary with good parameters for nice 
	publication figures.
	
	Resulting ``dict`` can be loaded via ``plt.rcParams.update()``.
	
	.. note:: Use this if you want to include LaTeX symbols in the figure.
	
	Keyword Args:
		fontSize (int): Font size.
		labelFontSize (int): Font size.
		tickFontSize (int): Font size.
		legendFontSize (int): Font size.
		
	Returns:
		dict: Parameter dictionary.
	
	"""
	
	
	params = {'backend': 'ps',
		'font.size': fontSize,
		'axes.labelsize': labelFontSize,
		'legend.fontsize': legendFontSize,
		'xtick.labelsize': tickFontSize,
		'ytick.labelsize': tickFontSize,
		'text.usetex': True,
		'font.family': 'sans-serif',
		'font.sans-serif': 'Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Arial, Helvetica, Avant Garde, sans-serif',
		#'ytick.direction': 'out',
		'text.latex.preamble': [r'\usepackage{helvet}', r'\usepackage{sansmath}',r'\usepackage{siunitx}'] , #r'\sansmath',    
		'ytick.direction' : 'out',
		'xtick.direction' : 'out'
		}
		
	return params

def turnAxesForPub(ax,adjustFigSize=True,figWidthPt=180.4,figHeightPt=None,ptPerInches=72.27,fontSize=10,labelFontSize=10,tickFontSize=10,legendFontSize=10,is3D=False):
	
	"""Turns axes nice for publication.
	
	If ``adjustFigSize=True``, will also adjust the size the figure. 
	
	Args:
		ax (matplotlib.axes): A matplotlib axes.
		
	Keyword Args:
		adjustFigSize (bool): Adjust the size of the figure.
		figWidthPt (float): Width of the figure in pt.
		figHeightPt (float): Height of the figure in pt.
		ptPerInches (float): Resolution in pt/inches.
		fontSize (int): Font size.
		labelFontSize (int): Font size.
		tickFontSize (int): Font size.
		legendFontSize (int): Font size.
		is3D (bool): Flag indicating if 3D figure.
		
	Returns:
		matplotlib.axes: Modified matplotlib axes.
		
	"""

	params=getPubParms(fontSize=fontSize,labelFontSize=labelFontSize,tickFontSize=tickFontSize,legendFontSize=legendFontSize)
	plt.rcParams.update(params)
	
	if not is3D:
		ax=setPubAxis(ax)
	
	setPubFigSize(ax.get_figure(),figWidthPt=figWidthPt,figHeightPt=figHeightPt,ptPerInches=ptPerInches)
	
	ax=closerLabels(ax,padx=3,pady=1)
	
	return ax
		
def setPubAxis(ax):		
	
	"""Gets rid of top and right axis.
	
	Args:
		ax (matplotlib.axes): A matplotlib axes.
		
	Returns:
		matplotlib.axes: Modified matplotlib axes.
	
	"""
	
	ax.spines['top'].set_color('none')
	ax.spines['right'].set_color('none')
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	
	ax.spines["left"].axis.axes.tick_params(direction="outward") 
	ax.spines["bottom"].axis.axes.tick_params(direction="outward") 
	
	return ax
		
def setPubFigSize(fig,figWidthPt=180.4,figHeightPt=None,ptPerInches=72.27):
	
	"""Adjusts figure size/aspect.
	
	If ``figHeightPt`` is not given, will use golden ratio to
	compute it.
	
	Keyword Args:
		figWidthPt (float): Width of the figure in pt.
		figHeightPt (float): Height of the figure in pt.
		ptPerInches (float): Resolution in pt/inches.
	
	Returns:
		matplotlib.figure: Adjust figure.
	
	"""
	
	inchesPerPt = 1.0/ptPerInches
	goldenMean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
	figWidth = figWidthPt*inchesPerPt  # width in inches
	if figHeightPt==None:
		figHeight = figWidth*goldenMean      # height in inches
	else:
		figHeight = figHeightPt*inchesPerPt
	figSize =  [figWidth,figHeight]
	
	fig.set_size_inches(figSize[0],figSize[1])
	
	fig.subplots_adjust(bottom=0.25)
	fig.subplots_adjust(left=0.2)
	fig.subplots_adjust(top=0.9)
	
	return fig

def closerLabels(ax,padx=10,pady=10):
	
	"""Moves x/y labels closer to axis."""
	
	ax.xaxis.labelpad = padx
	ax.yaxis.labelpad = pady
	return ax

def getRandomColor():
	
	"""Returns triplet defining a random color.
	"""
	
	return np.random.rand(3,1)

def is3DAxes(ax):
	
	"""Returns if an axes is a 3D axis.
	
	Args:
		ax (matplotlib.axes): A matplotlib axes.
		
	Returns:
		bool: True if 3d axis.
	"""
	
	if hasattr(ax, 'get_zlim'): 
		return True
	else:
		return False

