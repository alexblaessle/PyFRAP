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

"""Statistics module for PyFRAP toolbox, mainly used to evaluate goodness of fit, but also providing functions
to assess overall measurement statistics.
"""

#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

#Numpy
import numpy as np
import scipy

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

def computeFitRsq(fit):
	
	"""Computes R-squared values for fit object.
	
	R-squared values contain:
		
		* Mean R-squared value over all ROIs included in fit, stored in ``fit.MeanRsq``.
		* Product of R-squared value over all ROIs included in fit, stored in ``fit.Rsq``.
		* R-squared value for each ROIs included in fit, stored in ``fit.RsqBuROI``.
		
	Args:
		fit (pyfrp.subclasses.pyfrp_fit.fit): Fit object.
		
	Returns:
		pyfrp.subclasses.pyfrp_fit.fit: Updated fit object.
	
	"""
	
	Rsqs=[]
	
	#Compute Rsqs for regions used for fitting	
	for i in range(len(fit.dataVecsFitted)):
		r=Rsq(fit.dataVecsFitted[i],fit.fittedVecs[i])
		fit.RsqByROI[fit.ROIsFitted[i].name]=r
		Rsqs.append(r)
	
	fit.MeanRsq=np.mean(Rsqs)	
	fit.Rsq=np.prod(Rsqs)
		
	return fit

def Rsq(data,x):
	
	r"""Computes R-squared values for fit series to data series.
	
	R-squared value is being computed as 
	
	.. math:: R^2 = 1 - \frac{\sum\limits_i (x_i - d_i)^2}{\sum\limits_i (d_i - \bar{d} )^2}
	
	Args:
		x (numpy.ndarray) Fit series.
		data (numpy.ndarray): Data series.
		
	Returns:
		float: R-squared value.
	
	"""
	
	#Make numpy array to avoid list substraction problem
	data=np.asarray(data)
	x=np.asarray(x)
	
	#Compute ssd between x and data
	ssd=computeSSD(data,x)
	
	#Compute mean data
	meanData=np.mean(data)
	
	#Get derivation of datapoints from mean
	SStot=sum((data-meanData)**2)
	
	#Calculate RSq value
	Rsq=1-ssd/SStot
	
	return Rsq

def computeSSD(data,x):
	
	r"""Computes sum of squared differences (SSD) of fit series to data series.
	
	The SSD is computed by
	
	.. math:: SSD = \sum\limits_i (x_i - d_i)^2
	
	Args:
		x (numpy.ndarray) Fit series.
		data (numpy.ndarray): Data series.
		
	Returns:
		float: SSD.
	"""
	
	return sum((data-x)**2)

def parameterStats(x):
	
	r"""Returns mean, standard deviation and standard error of array.
	
	Note that standard error is computed as 
	
	.. math:: \sigma_n = \frac{\sigma}{\sqrt{n}}
	
	where :math:`n` is the number of samples in ``x`` and :math: :math:`\sigma`
	is the standard deviation over ``x``.
	
	Args:
		x (numpy.ndarray): List of values.
	
	Returns:
		tuple: Tuple containing:
			
			* xMean (float): Mean of x.
			* xStd (float): Standard deviation of x.
			* xSterr (float): Standard error of x.
			
	"""
	
	return np.mean(x), np.std(x), np.std(x)/np.sqrt(len(x))
	
def overlapSubSampleSelect(d,n,k):
	
	r"""Takes subsamples of size ``n`` that overlap on both 
	sides by ``k`` points out of d.
	
	Algorithm collects snippets of d from :math:`j(n-k)` to
	:math:`(j+1)n-jk`, where :math:`j` is the counter for the
	subsamples.
	
	Args:
		d (numpy.ndarray): Data vector.
		n (int): Size of subsamples.
		k (int): Overlap.
		
	Returns:
		list: List of ``numpy.ndarray`` of overlapping subsamples.
		
	"""
	
	j=0
	ds=[]
	for i in range(len(d)):
		try:
			dnew=d[j*(n-k):(j+1)*n-j*k]
		except IndexError:
			dnew=d[j*(n-k):]
		ds.append(dnew)
		
		if (j+1)*n-j*k>len(d):
			break
		
		j=j+1
		
	return ds	

def overlapSubSampleError(d,n,k):
	
	r"""Computes error between overlapping subsamples.
	
	Error is calculated by
	
	.. math:: \left|\frac{\bar{d_i}+\epsilon}{\bar{d_j}+\epsilon}\right|,
	
	where :math:`i,j \in {1,..,\frac{N}{n-k}}` and :math:`\epsilon` is some 
	offset to avoid singularties.
	
	.. note:: The resulting error matrix is symmetric.
	
	Args:
		d (numpy.ndarray): Data vector.
		n (int): Size of subsamples.
		k (int): Overlap.
		
	Returns:
		numpy.ndarray: Error matrix.
	
	"""
	
	ds=overlapSubSampleSelect(d,n,k)
	
	dMean=[]
	for dnew in ds:
		dMean.append(np.mean(dnew))
	
	dError=np.zeros((len(dMean),len(dMean)))
	for i in range(len(dMean)):
		for j in range(len(dMean)):
			dError[i,j]=abs((dMean[i]+1E-10)/(dMean[j]+1E-10+1E-11)-1)
		
	return dError	
	
def selectDataByOverlapSubSample(d,n,k,thresh,debug=False):
	
	"""Selects data vector based on overlapping
	subsampling and simple thresholding.
	
	This algorithm combines local derivatives with global changes
	and filters both datasets that have large local changes as well
	as large global changes. However, taking means over subsamples
	prevents neglecting data sets that have short peaks over only 1
	or 2 time points, such as bubbles etc.
	
	Args:
		d (numpy.ndarray): Data vector.
		n (int): Size of subsamples.
		k (int): Overlap.
		thresh (float): Selecting threshold.
		
	Keyword Args:
		debug (bool): Print debugging messages.
		
	Returns:
		bool: True if data set should be neglected.
	
	"""

	dError=overlapSubSampleError(d,n,k)		
	
	idxs=np.where(dError>thresh)
	
	if debug:
		
		if len(dError[idxs].flatten())>0:
			for i in range(len(idxs[0])):
				print("Subsamples ", idxs[0][i],idxs[1][i], " generate error ", dError[idxs[0][i],idxs[1][i]])
			
	return len(dError[idxs].flatten())>0

def computeLogLikelihood(fit,ROIs=None,neg=True,sigma=1):
	
	r"""Computes log-likelihood of fit assuming
	normal distribution of data around fit.
	
	Generally we assume that the data is distributed normally around
	the fitted model, that is
	
	.. math:: (m-d) \sim \mathcal{N}(0,\sigma)
	
	Assuming that all measurements are independent, this results in 
	a likelihood function of
	
	.. math:: \prod\limits_{j=1}^{n} \frac{1}{\sigma \sqrt{2\pi}} \exp\left(-\frac{(m_j-d_j)^2}{2\sigma^2}\right) 
	
	resulting in a log-likelihood function of
	
	.. math:: - \frac{n}{2} \log (2\pi) - \frac{\sum\limits_{j=1}^{n} (m_j-d_j)^2 }{2\sigma^2}
	
	Function allows to compute both negative and positive log-likelihood by ``neg`` flag.
	
	If the standard deviation :math:`sigma` is unknown, but the log-likelihood is
	only needed for model comparison reasons over the same dataset, all `sigma` dependent
	terms can be ignored. Thus this function returns both SSD and log-likelihood.
	
	.. note:: If ``ROIs=None``, will use all ROIs defined in ``fit.ROIsFitted``. If ``ROIs`` are
	   specified, but not in ``fit.ROIsFitted``, will simply skip them.
	
	Args:
		fit (pyfrp.subclasses.pyfrp_fit.fit): Fit object.
		
	Keyword Args:	
		ROIs (list): List of ROIs to be considered for computation.
		neg (bool): Compute negative log-likehood.
		sigma (float): Standard deviation of normal distribution if known.
		
	Returns:
		tuple: Tuple containing:
		
			* sign (float): Sign of likelihood function (-1/1).
			* SSD (float): SSD of fit over all ROIs specified in ROIs.
			* logLL (float): Log-likehood function at ``sigma=sigma``
	
	"""
	
	#Recompute SSD
	ssds=[]
	for i,r in enumerate(fit.ROIsFitted):
		if ROIs!=None:
			if r not in ROIs:
				continue
		
		ssds.append(computeSSD(fit.dataVecsFitted[i],fit.fittedVecs[i]))
	
	SSD=sum(ssds)
	
	#Compute log likelihood
	n = len(fit.dataVecsFitted[0])
	sign = (2*int(neg)-1)
	const = n/2.*np.log(2*np.pi)
	
	logLL=sign * ( const + SSD/(2*sigma))
	
	return sign,SSD,logLL

def computeAIC(fit,ROIs=None,sigma=1,fromSSD=True):
	
	r"""Computes Akaike Information Criterion of fit.
	
	The AIC is defined by
	
	.. math:: AIC= 2k - 2\log(L),
	
	:math:`k` is the number of free parameters of the model used in 
	``fit`` and :math:`L` is the maximum likelihood function of the
	model, see also :py:func:`computeLogLikelihood`.
	
	Since generally the AIC is only used for model comparison and a normal
	distribution of the data around the model is assumed, one simply
	can use the SSD as a log-likelihood function. ``fromSSD`` controls if
	just the SSD is used. This is also useful, since :math:`sigma` is not 
	necessarily known if the objective function of the fit was not a
	maximum likelihood function.
	
	.. note:: If ``ROIs=None``, will use all ROIs defined in ``fit.ROIsFitted``. If ``ROIs`` are
	   specified, but not in ``fit.ROIsFitted``, will simply skip them.
	   
	If the AIC or the AICc should be used can be determined using :py:func:`useAIC`.
	
	Args:
		fit (pyfrp.subclasses.pyfrp_fit.fit): Fit object.
		
	Keyword Args:	
		ROIs (list): List of ROIs to be considered for computation.
		sigma (float): Standard deviation of normal distribution if known.
		fromSSD (bool): Simply use SSD as maximum likelihood.
		
	Returns:
		float: AIC of fit.
	
	"""
		
	AIC=2*(fit.getNParmsFitted())+2*computeLogLikelihood(fit,ROIs=ROIs)[2-int(fromSSD)]
	
	return AIC

def computeCorrAIC(fit,ROIs=None,sigma=1,fromSSD=True):
	
	r"""Computes corrected Akaike Information Criterion of fit.
	
	The AICc is defined by
	
	.. math:: AICc= AIC + \frac{2k(k+1)}{n-k-1},
	
	where :math:`n` is the number of datapoints used for the fit and 
	:math:`k` is the number of free parameters. For the computation of the
	AIC, please refer to the documentation of :py:func:`computeAIC`.
	
	.. note:: If ``ROIs=None``, will use all ROIs defined in ``fit.ROIsFitted``. If ``ROIs`` are
	   specified, but not in ``fit.ROIsFitted``, will simply skip them.
	   
	If the AIC or the AICc should be used can be determined using :py:func:`useAIC`.
	
	Args:
		fit (pyfrp.subclasses.pyfrp_fit.fit): Fit object.
		
	Keyword Args:	
		ROIs (list): List of ROIs to be considered for computation.
		sigma (float): Standard deviation of normal distribution if known.
		fromSSD (bool): Simply use SSD as maximum likelihood.
		
	Returns:
		float: AICc of fit.
	
	"""
	
	n = len(fit.dataVecsFitted[0])
	AIC=computeAIC(fit,ROIs=ROIs,sigma=sigma,fromSSD=fromSSD)
	k=fit.getNParmsFitted()
	
	corrAIC=AIC + 2*k*(k+1)/(n-k-1)
	
	return corrAIC
	
def useAIC(n,k):

	r"""Returns if corrected or not corrected version of the 
	AIC should be used.
	
	Rule of thumb is that AIC should be used if
	
	.. math:: \frac{n}{k}>40,
	
	where :math:`n` is the number of datapoints used for the fit and 
	:math:`k` is the number of free parameters.
	
	Returns:
		bool: ``True``, use AIC, ``False`` use AICc.
	
	"""
	
	return n/k>40

def compareFitsByAIC(fits,ROIs=None,sigma=1,fromSSD=True,thresh=None):
	
	r"""Computes AIC and Akaike weights for all fits in a list
	and returning best models given certain criteria.
	
	For the computation of the AIC see :py:func:`computeAIC` and 
	the computation of the Akaike weights :py:func:`computeAkaikeWeights`.
	
	If threshold ``thresh=None``, then will select model with maximum
	Akaike weight as best model, that is, the model with the highest likelihood
	of being the best model.
	If ``thresh`` is given, will return list of acceptable models based on
	
	.. math:: w_{\mathrm{max}}-w_i<thresh,
	
	that is, all models that are within a given range of probability of
	the most likely one.
	
	Args:
		fits (list): List of fit objects.
		
	Keyword Args:	
		ROIs (list): List of ROIs to be considered for computation.
		sigma (float): Standard deviation of normal distribution if known.
		fromSSD (bool): Simply use SSD as maximum likelihood.
		thresh (float): Probability range for model selection.
	
	Returns:
		tuple: Tuple containing:
		
			* AICs (list): List of AIC values of the respective fits.
			* deltaAICs (numpy.ndarray): List of Akaike difference values of the respective fits.
			* weights (numpy.ndarray): List of Akaike difference weights of the respective fits.
			* acc (list): List of acceptable fits by model selection.
			* ks (list): List of number of parameters fitted of the respective fits.
			* ns (list): List of number of datapoints fitted of the respective fits.
			
	"""
	
	AICs=[]
	ks=[]
	ns=[]
	for fit in fits:
		AICs.append(computeAIC(fit,ROIs=ROIs,sigma=1,fromSSD=fromSSD))
		ks.append(fit.getNParmsFitted())
		ns.append(len(fit.dataVecsFitted[0]))
		
	deltaAIC=np.asarray(AICs)-min(AICs)	
	
	weights = computeAkaikeWeights(AICs)
	
	acc=[]
	if thresh==None:
		
		idx=list(weights).index(max(weights))
		
		acc=[fits[idx]]
	
	else:
		relWeights=abs(weights-max(weights))
		for i in np.where(list(relWeights<thresh))[0]:
			acc.append(fits[i])
	
	
	return AICs, deltaAIC,weights, acc,ks,ns

def compareFitsByCorrAIC(fits,ROIs=None,sigma=1,fromSSD=True,thresh=None):
	
	r"""Computes AICc and Akaike weights for all fits in a list
	and returning best models given certain criteria.
	
	For the computation of the corrected AIC see :py:func:`computeCorrAIC` and 
	the computation of the Akaike weights :py:func:`computeAkaikeWeights`.
	
	If threshold ``thresh=None``, then will select model with maximum
	Akaike weight as best model, that is, the model with the highest likelihood
	of being the best model.
	If ``thresh`` is given, will return list of acceptable models based on
	
	.. math:: w_{\mathrm{max}}-w_i<thresh,
	
	that is, all models that are within a given range of probability of
	the most likely one.
	
	Args:
		fits (list): List of fit objects.
		
	Keyword Args:	
		ROIs (list): List of ROIs to be considered for computation.
		sigma (float): Standard deviation of normal distribution if known.
		fromSSD (bool): Simply use SSD as maximum likelihood.
		thresh (float): Probability range for model selection.
	
	Returns:
		tuple: Tuple containing:
		
			* AICs (list): List of AICc values of the respective fits.
			* deltaAICs (numpy.ndarray): List of Akaike difference values of the respective fits.
			* weights (numpy.ndarray): List of Akaike difference weights of the respective fits.
			* acc (list): List of acceptable fits by model selection.
			* ks (list): List of number of parameters fitted of the respective fits.
			* ns (list): List of number of datapoints fitted of the respective fits.
	
	"""
	
	AICs=[]
	ks=[]
	ns=[]
	for fit in fits:
		AICs.append(computeCorrAIC(fit,ROIs=ROIs,sigma=1,fromSSD=fromSSD))
		ks.append(fit.getNParmsFitted())
		ns.append(len(fit.dataVecsFitted[0]))
		
	deltaAIC=np.asarray(AICs)-min(AICs)	
	
	weights = computeAkaikeWeights(AICs)
	
	acc=[]
	if thresh==None:
		
		idx=list(weights).index(max(weights))
		acc=[fits[idx]]
	
	else:
		relWeights=abs(weights-max(weights))
		
		for i in np.where(list(relWeights<thresh))[0]:
			acc.append(fits[i])
	
	return AICs, deltaAIC,weights, acc,ks,ns

	
def computeAkaikeWeights(AICs):
	
	r"""Computes Akaike weights for a list of 
	AIC values.
	
	Akaike weightsare given by:
	
	.. math:: w_i = \frac{\exp\left(\frac{-\Delta_{i}}{2}\right)}{\sum\limits_{i=1}^{N} \exp\left(\frac{- \Delta_{i}}{2}\right)},
	
	where :math:`\Delta_{i}` is the Akaike difference of model :math:`i` computed by
	
	.. math:: \Delta_{i} = AIC_i - AIC_{\mathrm{min}},
	
	and :math:`N` is the total number of models considered. Note that
	
	.. math:: \sum\limits_{i=1}^{N} \exp\left(\frac{- \Delta_{i}}{2}\right) = 1.
	
	Args:
		AICs (list): List of AIC values.
		
	Returns:
		np.ndarray: Corresponding Akaike weights.
	
	
	"""
	
	deltaAIC=np.asarray(AICs)-min(AICs)	
	
	weights=np.exp(-deltaAIC/2.)/sum(np.exp(-deltaAIC/2.))
	
	
	return weights

def wilcoxonTest(x,y,zero_method='wilcox', correction=False,printOut=True):
	
	"""Performs wilcoxon test between two groups.
	
	See also https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.wilcoxon.html
	
	Args:
		x (list): First group.
		y (list): Second group.
		
	Keyword Args:
		zero_method (str): Treatment of zeros.
		correction (bool): Apply continuety correction.
		printOut (bool): Print out results.
	
	Returns:
		tuple: Tuple containing:
		
			* stat (float): Sum of ranks.
			* val (float): p-Value.
			
	"""
	
	stat,pval=scipy.stats.wilcoxon(x, y=y, zero_method=zero_method, correction=correction)
	
	if printOut:
		print("Results of Wilcoxon-Test:")
		print("p-Value: ", pval)
		print("Wilcoxon-Statistics:", stat)
	
	return stat,pval

def mannWhitneyUTest(x,y,printOut=True):
	
	"""Performs Mann-Whitney-U test between two groups.
	
	See also https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.mannwhitneyu.html#scipy.stats.mannwhitneyu
	
	Args:
		x (list): First group.
		y (list): Second group.
		
	Keyword Args:
		printOut (bool): Print out results.
	
	Returns:
		tuple: Tuple containing:
		
			* stat (float): Mann-Whitney statistics.
			* val (float): p-Value.
			
	"""
	
	stat,pval=scipy.stats.mannwhitneyu(x, y)
	
	if printOut:
		print("Results of Mann-Whitney-U-Test:")
		print("p-Value: ", pval)
		print("U-Statistics:", stat)
	
	return stat,pval

def tTestStandard(x,y,printOut=True):
	
	"""Performs standard t-test between two groups.
	
	See also https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.ttest_ind.html#scipy.stats.ttest_ind
	
	Args:
		x (list): First group.
		y (list): Second group.
		
	Keyword Args:
		printOut (bool): Print out results.
	
	Returns:
		tuple: Tuple containing:
		
			* stat (float): Test statistics.
			* val (float): p-Value.
			
	"""
	
	stat,pval=scipy.stats.ttest_ind(x, y, equal_var=True)
	
	if printOut:
		print("Results of Standard t-Test:")
		print("p-Value: ", pval)
		print("Statistics:", stat)
		
	return stat,pval
	
def tTestWelch(x,y,printOut=True):
	
	"""Performs Welch's t-test between two groups.
	
	See also https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.ttest_ind.html#scipy.stats.ttest_ind
	
	Args:
		x (list): First group.
		y (list): Second group.
		
	Keyword Args:
		printOut (bool): Print out results.
	
	Returns:
		tuple: Tuple containing:
		
			* stat (float): Test statistics.
			* val (float): p-Value.
			
	"""
	
	stat,pval=scipy.stats.ttest_ind(x, y, equal_var=False)
	
	if printOut:
		print("Results of Welch's t-Test:")
		print("p-Value: ", pval)
		print("Statistics:", stat)	
	
	return stat,pval
	
def shapiroTest(x,printOut=True):
	
	"""Performs shapiro test to test group for normality.
	
	See also https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.shapiro.html#scipy.stats.shapiro
	
	Args:
		x (list): Data.
		
	Keyword Args:
		printOut (bool): Print out results.
	
	Returns:
		tuple: Tuple containing:
		
			* stat (float): Test statistic.
			* val (float): p-Value.
			
	"""

	stat,pval=scipy.stats.sharipo(x)
	
	if printOut:
		print("Results of Sharipo-Test:")
		print("p-Value: ", pval)
		print("Statistics:", stat)
	
	return stat,pval

	
	
	
	
	
	
	