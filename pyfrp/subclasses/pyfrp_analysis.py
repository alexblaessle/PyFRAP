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

"""Essential PyFRAP module containing :py:class:`analysis` class. 
"""

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#Numpy/Scipy
from numpy import *

#JSON Compression
import pickle

#PyFRAP Modules
from pyfrp.modules import pyfrp_img_module 
from pyfrp.modules import pyfrp_misc_module
from pyfrp.modules.pyfrp_term_module import *

#matplotlib
import matplotlib.pyplot as plt

#===========================================================================================================================================================================
#Module Classes
#===========================================================================================================================================================================


class analysis:
	
	"""PyFRAP analysis class storing information about analysis options 
	and some analysis results.
	
	Analysis options are:
	
		* ``gaussian``: Apply gaussian filter to images. Default kernel size is ``gaussianSigma=2``.
		* ``median``: Apply gaussian filter to images. Default kernel size is ``medianRadius=5``.
		* ``flatten``: Apply flattening mask.
		* ``norm``: Norm  by pre image.
		* ``bkgd``: Substract background.
		* ``quad``: Perform reduction to first quadrant by flipping.
		* ``flipBeforeProcess``: Flip into quadrant before other processing options are applied.
	
	Analysis options are stored in ``process`` dictionary. If analysis finds option in ``process.keys``, it will
	perform option. Analysis options can be turned on/off using the respective functions, such as 
	
		* :py:func:`pyfrp.subclasses.pyfrp_analysis.medianOn`
		* :py:func:`pyfrp.subclasses.pyfrp_analysis.flattenOn`
		* etc.
	
	Processing parameters are stored in ``process.values``.
	
	The default processing options are ``process={}``, meaning that no image modification is applied before 
	concentration readout, see also :py:func:`genDefaultProcess`.
	
	.. warning:: Quadrant reduction is still experimental.
	
	Three other important attributes are:
	
		* ``dataOffset``: The offset of the data that is for example used for norming, see also :py:func:`getOptimalOffset`.
		* ``addRimImg``: Flag that controls if rim concentrations are added to ROI concentratrtion profiles, see also
		  :py:func:`setAddRimImg`.
		* ``concRim``: The rim concentration of the first post-bleaching image used later by the simulation for nodes that 
		  are outside of original image boundaries.
	
	.. note:: ``addRimImg=True`` by default. This is generally good, since the simulation value in ROIs is getting evaluated over
	   over mesh nodes both inside the actual image and outside of it. 
	
	Args:
		embryo (pyfrp.subclasses.pyfrp_embryo.embryo): PyFRAP embryo instance.
	
	"""
	
	#Creates new embryo object
	def __init__(self,embryo):
		
		self.embryo=embryo
		
		#Rim handling
		self.addRimImg=True
		self.concRim=None
		
		#Norming Data
		self.fnPreimage=''
		self.fnFlatten=''
		self.fnBkgd=''
		self.flatteningMask=None
		self.bkgdMask=None
		self.preMask=None
		self.nPre=2
		self.nFlatten=2
		self.nBkgd=2
		
		#Data offset
		self.dataOffset=1.
			
		#Processing options
		self.process={}
		self.gaussianSigma=2
		self.medianRadius=5
		
	def run(self,signal=None,embCount=None,debug=False,debugAll=False,showProgress=True):
		
		"""Runs analysis by passing analysis object to :py:func:`pyfrp.modules.pyfrp_img_module.analyzeDataset`.
		
		Will first check if ROI indices are computed for all ROIs and if necessary compute them before starting 
		data analysis.
		
		Keyword Args:
			signal (PyQt4.QtCore.pyqtSignal): PyQT signal to send progress to GUI.
			embCount (int): Counter of counter process if multiple datasets are analyzed. 
			debug (bool): Print final debugging messages and show debugging plots.
			debugAll (bool): Print debugging messages and show debugging plots of each step.
			showProgress (bool): Print out progress.
		
		Returns:
			pyfrp.subclasses.pyfrp_analysis.analysis: Updated analysis instance.
		
		"""
		
		if 'norm' in self.process and 'flatten' in  self.process:
			printWarning("Both norm and flatten have been selected for data analysis. This is not advisable.")
		
		if not self.embryo.checkROIIdxs()[0]:
			self.embryo.computeROIIdxs()
			
		self=pyfrp_img_module.analyzeDataset(self,signal=signal,embCount=embCount,debug=debug,debugAll=debugAll,showProgress=showProgress)
		return self
	
	def setGaussianSigma(self,s):
		
		"""Sets size of gaussian kernel and updates its value
		in ``process`` dictionary if gaussian filter is turned on.
		
		See also http://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.gaussian_filter.
		
		Args:
			s (float): New sigma.
			
		"""
		
		self.gaussianSigma=s
		self.updateProcess()
		return self.gaussianSigma
	
	def getGaussianSigma(self):
		
		"""Returns size of gaussian kernel.
		
		See also http://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.gaussian_filter.
		
		Returns:
			float: Gaussian sigma.
			
		"""
		
		return self.gaussianSigma
	
	def setMedianRadius(self,s):
		
		"""Sets size of median kernel and updates its value
		in ``process`` dictionary if median filter is turned on.
		
		See also http://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.median and 
		http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.ndimage.filters.median_filter.html.
		
		Args:
			s (float): New radius.
			
		"""
		
		self.medianRadius=s
		self.updateProcess()
		return self.medianRadius
	
	def getMedianRadius(self):
		
		"""Returns size of median kernel.
		
		See also http://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.median and 
		http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.ndimage.filters.median_filter.html.
		
		Returns:
			float: New radius.
			
		"""
		
		return self.medianRadius
	
	def setGaussian(self,b):
		
		"""Turns on/off gaussian filter for analysis.
		
		.. note:: Will use ``gaussianSigma`` as kernel size. Can be changed via 
		   :py:func:`setGaussianSigma`.
		
		Args:
			b (bool): ``True`` if gaussian should be turned on, ``False`` else.
			
		Returns:
			dict: Updated process dictionary.
		
		"""
		
		return self.parm2Process(b,'gaussian',self.gaussianSigma)
		
	def setNorm(self,b):
		
		"""Turns on/off norming by preimage for analysis.
		
		.. note:: Will use ``preMask`` for norming. ``preMask`` is updated
		   via :py:func:`computePreMask` and then automatically updated in ``process`` dictionary.
		
		Args:
			b (bool): ``True`` if norming should be turned on, ``False`` else.
			
		Returns:
			dict: Updated process dictionary.
		
		"""
		
		return self.parm2Process(b,'norm',self.fnPreimage)
	
	def setFlipBeforeProcess(self,b):
		
		"""Turns on/off if image should be flipped into quadrant before or after 
		performing all other image processing for analysis.
		
		.. warning:: Quadrant reduction is still experimental.
		
		Args:
			b (bool): ``True`` if image should be flipped before, ``False`` else.
			
		Returns:
			dict: Updated process dictionary.
		
		"""
		
		return self.parm2Process(b,'flipBeforeProcess',True)
		
	def setMedian(self,b):
		
		"""Turns on/off median filter for analysis.
		
		.. note:: Will use ``medianRadius`` as kernel size. Can be changed via 
		   :py:func:`setMedianRadius`.
		
		Args:
			b (bool): ``True`` if median should be turned on, ``False`` else.
			
		Returns:
			dict: Updated process dictionary.
		
		"""
		
		return self.parm2Process(b,'median',self.medianRadius)
	
	def setQuad(self,b):
		
		"""Turns on/off if image should be flipped into first quadrant for analysis.
		
		.. warning:: Quadrant reduction is still experimental.
		
		Args:
			b (bool): ``True`` if quadrant reduction should be turned on, ``False`` else.
			
		Returns:
			dict: Updated process dictionary.
		
		"""
		
		return self.parm2Process(b,'quad',True)
	
	def setFlatten(self,b):
		
		"""Turns on/off flattening for analysis.
		
		.. note:: Will use ``flatteningMask`` for flattening. ``flatteningMask`` is updated
		   via :py:func:`computeFlatteningMask` and then automatically updated in ``process`` dictionary.
		
		Args:
			b (bool): ``True`` if flattening should be turned on, ``False`` else.
			
		Returns:
			dict: Updated process dictionary.
		
		"""
		
		return self.parm2Process(b,'flatten',self.fnFlatten)
	
	def setBkgd(self,b):
		
		"""Turns on/off background substraction for analysis.
		
		.. note:: Will use ``bkgdMask`` for flattening. ``bkgdMask`` is updated
		   via :py:func:`computeBkgdMask` and then automatically updated in ``process`` dictionary.
		
		Args:
			b (bool): ``True`` if background substraction should be turned on, ``False`` else.
			
		Returns:
			dict: Updated process dictionary.
		
		"""
		
		return self.parm2Process(b,'bkgd',self.fnBkgd)
	
	def setFnFlatten(self,fn):
		
		"""Sets path to flattening dataset.
		
		Args:
			fn (str): Path to flattening dataset.
		
		"""
		
		self.fnFlatten=fn
		self.updateProcess()
		return self.fnFlatten
	
	def getFnFlatten(self):
		
		"""Returns path to flattening dataset.
		
		Returns:
			str: Path to flattening dataset.
		"""
		
		return self.fnFlatten
	
	def setFnBkgd(self,fn):
		
		"""Sets path to background dataset.
		
		Args:
			fn (str): Path to background dataset.
		
		"""
		
		self.fnBkgd=fn
		self.updateProcess()
		return self.fnBkgd
	
	def getFnBkgd(self):
			
		"""Returns path to background dataset.
		
		Returns:
			str: Path to background dataset.
		"""
		
		return self.fnBkgd
	
	def setFnPre(self,fn):
		
		"""Sets path to preimage dataset.
		
		Args:
			fn (str): Path to preimage dataset.
		
		"""
		
		self.fnPreimage=fn
		self.updateProcess()
		return self.fnPreimage
	
	def getFnPre(self):
		
		"""Returns path to norming dataset.
		
		Returns:
			str: Path to norming dataset.
		"""
		
		return self.fnPreimage
	
	def setDataOffset(self,s):
		
		"""Sets dataoffset used for norming.
		
		Args:
			s (float): New offset.

		"""
		
		self.dataOffset=s 
		return self.dataOffset
	
	def getDataOffset(self):
		
		"""Returns dataoffset used for norming.
		
		Returns:
			float: Current offset.

		"""
		
		return self.dataOffset
	
	def setAddRimImg(self,s):
		
		"""Sets the addRimImg flag.
		
		The addRim flag controls if the rim concentration is added to the concentration 
		of each ROI timeseries depending on how many imaginary pixels they have outside
		of the actual image.
		
		Args:
			s (bool): Flag value.
			
			
		"""
		
		self.addRimImg=s 
		return self.addRimImg
	
	def getAddRimImg(self):
		
		"""Returns the addRimImg flag.
		
		See also :py:func:`setAddRimImg`.
		
		Returns:
			bool: Flag value.
			
		"""
		
		return self.addRimImg
	
	def setConcRim(self,s):
		
		"""Sets rim concentration.
		
		Args:
			s (float): New rim concentration.
		"""
		
		self.concRim=s 
		return self.concRim
	
	def getConcRim(self):
		
		"""Returns rim concentration.
		
		Returns:
			float: Current rim concentration.
		"""
		
		return self.concRim
	
	def setProcess(self,s):
		
		"""Sets process dictionary.
		
		Args:
			s (dict): New process dictionary.
		"""
		
		self.process=s 
		return self.process
	
	def getProcess(self):
		
		"""Returns process dictionary."""
		
		return self.process
	
	def setNPre(self,n):
		
		"""Sets the number of images used for the computation
		for the mean norming image.
		
		Args:
			n (int): Number of images used.
			
		"""	
		
		self.nPre=n 
		return self.nPre
	
	def getNPre(self):
		
		"""Returns the number of images used for the computation
		for the mean norming image.
		
		Returns:
			int: Number of images used.
			
		"""	
		
		return self.nPre
	
	def setNFlatten(self,n):
		
		"""Sets the number of images used for the computation
		for the mean flattening image.
		
		Args:
			n (int): Number of images used.
			
		"""	
		
		self.nFlatten=n 
		return self.nFlatten
	
	def getNFlatten(self):
		
		"""Returns the number of images used for the computation
		for the mean flattening image.
		
		Returns:
			int: Number of images used.
			
		"""	
		
		return self.nFlatten
	
	def setNBkgd(self,n):
		
		"""Sets the number of images used for the computation
		for the mean background image.
		
		Args:
			n (int): Number of images used.
			
		"""	
		
		self.nBkgd=n 
		return self.nBkgd
	
	def getNBkgd(self):
		
		"""Returns the number of images used for the computation
		for the mean background image.
		
		Returns:
			int: Number of images used.
			
		"""	
		
		
		return self.nBkgd
	
	def normOn(self):
		
		"""Returns current state of this option.
		
		Returns:
			bool: ``True`` if switched on, ``False`` else.
		
		"""
		
		return 'norm' in list(self.process.keys())
	
	def bkgdOn(self):
		
		"""Returns current state of this option.
		
		Returns:
			bool: ``True`` if switched on, ``False`` else.
		
		"""
		
		return 'bkgd' in list(self.process.keys())
	
	def flattenOn(self):
		
		"""Returns current state of this option.
		
		Returns:
			bool: ``True`` if switched on, ``False`` else.
		
		"""
		
		return 'flatten' in list(self.process.keys())
	
	def gaussianOn(self):
		
		"""Returns current state of this option.
		
		Returns:
			bool: ``True`` if switched on, ``False`` else.
		
		"""
		
		return 'gaussian' in list(self.process.keys())
	
	def medianOn(self):
		
		"""Returns current state of this option.
		
		Returns:
			bool: ``True`` if switched on, ``False`` else.
		
		"""
		
		return 'median' in list(self.process.keys())
	
	def quadOn(self):
		
		"""Returns current state of this option.
		
		Returns:
			bool: ``True`` if switched on, ``False`` else.
		
		"""
		
		return 'quad' in list(self.process.keys())
	
	def flipBeforeProcessOn(self):
		
		"""Returns current state of this option.
		
		Returns:
			bool: ``True`` if switched on, ``False`` else.
		
		"""
		
		return 'flipBeforeProcess' in list(self.process.keys())
	
	def updateProcess(self):
		
		"""Updates all values in process dictionary with
		the ones saved in attributes of analysis object."""
		
		if 'flatten' in list(self.process.keys()):
			self.process['flatten']=self.fnFlatten
		if 'norm' in list(self.process.keys()):
			self.process['norm']=self.fnPreimage
		if 'bkgd' in list(self.process.keys()):
			self.process['bkgd']=self.fnBkgd	
		if 'gaussian' in list(self.process.keys()):
			self.process['gaussian']=self.gaussianSigma
		if 'median' in list(self.process.keys()):
			self.process['median']=self.medianRadius	
			
		return
	
	def printProcess(self):
		
		"""Prints out current process options in a nicely formatted way."""
		
		printDict(self.process)
		return
	
	def printAllAttr(self):
		
		"""Prints out all attributes of analysis object.""" 
		
		print("Analysis of embryo ", self.embryo.name, " Details.")
		printAllObjAttr(self)
	
	def parm2Process(self,b,key,val):
		
		"""Adds/Removes a new option to ``process`` dictionary.
		
		Args:
			b (bool): Flag if process should be added or removed.
			key (str): Key of option to be added.
			val (any): Value of dictionary entry.
		
		Returns:
			dict: Updated ``process`` dictionary.
		
		"""
		
		if b:
			self.process[key]=val
		else:	
			try:
				self.process.pop(key)
			except KeyError:
				pass
		return self.process	
		
	def genDefaultProcess(self):
		
		"""Sets ``process`` dictionary to default options.
		
		Default options are:
		
			* ``gaussian=False``
			* ``median=False``
			* ``quad=False``
			* ``flatten=False``
			* ``bkgd=False``
			* ``norm=False``
			* ``flipBeforeProcess=True``
		
		Returns:
			dict: Updated ``process`` dictionary.
		
		"""
		
		self.setGaussian(False)
		self.setFlipBeforeProcess(True)
		self.setQuad(False)
		self.setNorm(False)
		self.setMedian(False)
		self.setFlatten(False)
		self.setBkgd(False)
		return self.process
	
	def removeProcessStep(self,dic,step):
		
		"""Removes process step from dictionary.
		
		Args:
			dic (dict): A dictionary.
			step (str): Key of step to be removed.
			
		Returns:
			dict: Updated dictionary.
		
		"""
		
		try:
			dic.pop(step)
		except KeyError:
			pass
		return dic
		
	def computeFlatteningMask(self,applyProcess=True):
		
		"""Computes flattening mask.
		
		Takes first ``nFlatten`` images in ``fnFlatten`` and computes mean image of these 
		images. Then, if ``applyProcess`` is selected, applies the selected process options 
		defined in ``process`` dictionary to it.
		
		.. note:: Will not apply process options ``norm``, ``flatten`` and ``bkgd`` to mean 
		   flattening image.
		   
		Keyword Args:
			applyProcess (bool): Apply processing options to flattening mask.
			
		Returns:
			numpy.ndarray: Flattening mask.
		
		"""
		
		fileList=pyfrp_misc_module.getSortedFileList(self.fnFlatten,self.embryo.dataFT)
		fileList=fileList[:self.nFlatten]
		meanImg=pyfrp_img_module.computeMeanImg(self.fnFlatten,fileList,self.embryo.dataEnc)
		
		if applyProcess:
		
			processDic=dict(self.process)
			processDic=self.removeProcessStep(processDic,'norm')
			processDic=self.removeProcessStep(processDic,'flatten')
			processDic=self.removeProcessStep(processDic,'bkgd')
			
			###NOTE: Also remove bkgd??? Having bkgd in there could lead to singularities with flattening norming? 
			meanImg=pyfrp_img_module.processImg(meanImg,processDic,None,None,None,dataOffset=self.dataOffset)
		
		self.flatteningMask=pyfrp_img_module.computeFlatMask(meanImg,self.dataOffset)
	
		return self.flatteningMask
	
	def computeBkgdMask(self,flatteningMask,applyProcess=True,applyFlatten=False):
		
		"""Computes background mask.
		
		Takes first ``nBkgd`` images in ``fnBkgd`` and computes mean image of these 
		images. Then, if ``applyProcess`` is selected, applies the selected process options 
		defined in ``process`` dictionary to it.
		
		.. note:: Will not apply process options ``norm`` and ``bkgd`` to mean 
		   background image.
		
		Args:
			flatteningMask (numpy.ndarray): Flattening mask.
		
		Keyword Args:
			applyProcess (bool): Apply processing options to background mask.
			applyFlatten (bool): Apply flattening to background mask.
			
		Returns:
			numpy.ndarray: Background mask.
		
		"""
		
		fileList=pyfrp_misc_module.getSortedFileList(self.fnBkgd,self.embryo.dataFT)
		fileList=fileList[:self.nBkgd]
		meanImg=pyfrp_img_module.computeMeanImg(self.fnBkgd,fileList,self.embryo.dataEnc)
		
		if applyProcess:
			
			processDic=dict(self.process)
			processDic=self.removeProcessStep(processDic,'norm')
			processDic=self.removeProcessStep(processDic,'bkgd')
			if not applyFlatten:
				processDic=self.removeProcessStep(processDic,'flatten')
			
			meanImg=pyfrp_img_module.processImg(meanImg,processDic,flatteningMask,None,None,dataOffset=self.dataOffset)
		
		self.bkgdMask=meanImg
		
		return self.bkgdMask
	
	def computePreMask(self,flatteningMask,bkgdMask,applyProcess=True):
		
		"""Computes norming mask.
		
		Takes first ``nPre`` images in ``fnPreimage`` and computes mean image of these 
		images. Then, if ``applyProcess`` is selected, applies the selected process options 
		defined in ``process`` dictionary to it.
		
		.. note:: Will not apply process option ``norm`` to mean 
		   background image.
		
		Args:
			flatteningMask (numpy.ndarray): Flattening mask.
			bkgdMask (numpy.ndarray): Background mask.
			
		Keyword Args:
			applyProcess (bool): Apply processing options to background mask.
			
		Returns:
			numpy.ndarray: Norming mask.
		
		"""
		
		fileList=pyfrp_misc_module.getSortedFileList(self.fnPreimage,self.embryo.dataFT)
		fileList=fileList[:self.nPre]
		meanImg=pyfrp_img_module.computeMeanImg(self.fnPreimage,fileList,self.embryo.dataEnc)
		
		if applyProcess:
		
			processDic=dict(self.process)
			processDic=self.removeProcessStep(processDic,'norm')
			
			#removeProcesses=['norm','median']
			
			#processDic=self.removeProcessStep(processDic,'median')
			
			imgPre=pyfrp_img_module.processImg(meanImg,processDic,flatteningMask,bkgdMask,None,dataOffset=self.dataOffset)
		
		self.preMask=imgPre
		
		return imgPre
 
	def getOptimalOffset(self,debug=False):
		
		r"""Computes optimal dataoffset for data analysis.
		
		Finds minimal non-zero offset for main dataset, preimage dataset
		flattening dataset and background dataset, if available. The Idea is that
		one does not want to have negative pixels, so substraction of a fixed value
		from an image should always lead to positive pixel values. Thus the offset is 
		computed by
		
		.. math:: offset = max\{ o_{\mathrm{min},\mathrm{data}},o_{\mathrm{min},\mathrm{flatten}},o_{\mathrm{min},\mathrm{pre}},o_{\mathrm{min},\mathrm{bkgd}}\}
		
		where :math:`o_{\mathrm{min},d}` is the minimum pixel values of all images in dataset 
		:math:`d`.
		
		Keyword Args:
			debug (bool): Print debugging messages.
			
		Returns:
			float: Optimal offset.
			
		"""
		
		self.dataOffset=pyfrp_img_module.findMinOffset(self.embryo.fnDatafolder,self.embryo.fileList,self.embryo.dataEnc,oldOffset=self.dataOffset,defaultAdd=1.,debug=debug)
		
		if self.fnPreimage!=None:
			fileList=pyfrp_misc_module.getSortedFileList(self.fnPreimage,self.embryo.dataFT)[:self.nPre]
			self.dataOffset=pyfrp_img_module.findMinOffset(self.fnPreimage,fileList,self.embryo.dataEnc,oldOffset=self.dataOffset,defaultAdd=1.,debug=debug)
		
		if self.fnFlatten!=None:
			fileList=pyfrp_misc_module.getSortedFileList(self.fnFlatten,self.embryo.dataFT)[:self.nFlatten]
			self.dataOffset=pyfrp_img_module.findMinOffset(self.fnFlatten,fileList,self.embryo.dataEnc,oldOffset=self.dataOffset,defaultAdd=1.,debug=debug)
			
		if self.fnBkgd!=None:
			fileList=pyfrp_misc_module.getSortedFileList(self.fnBkgd,self.embryo.dataFT)[:self.nBkgd]
			self.dataOffset=pyfrp_img_module.findMinOffset(self.fnBkgd,fileList,self.embryo.dataEnc,oldOffset=self.dataOffset,defaultAdd=1.,debug=debug)
			
		return self.dataOffset
		