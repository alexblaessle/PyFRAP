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

"""Miscellaneous module for PyFRAP toolbox. Contains functions handling:

	* Filepath management (Windowns/Unix conversion etc.)
	* String searching/replacing
	* List matching/comparison
	* Dictionary to object and back conversion/extraction.
	* Settings path management.
	* Embryo wizard.
	* etc.
	
"""

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#PyFRAP modules
from .pyfrp_term_module import *
from pyfrp.modules import pyfrp_img_module

#PyFRAP subclasses


#Numpy
import numpy as np

#Misc
import csv
import time
import os
import inspect
import platform
import shutil
import glob
from tempfile import mkstemp

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

def getSortedFileList(fnFolder,fType):
	
	"""Gets sorted file list from folder for files of type fType
	
	Args:
		fnFolder (str): Folder path.
		fType (str): File type.
		
	Returns:
		list: list of filenames.
		
	"""
	
	#Getting files in datafolder
	try:
		fileList=os.listdir(fnFolder)
	except OSError:
		printWarning("Folder " + fnFolder + " does not exist.")
		return []
		
	fileListNew=[]
	
	#Going through all filenames 
	for i in range(len(fileList)):
		#Check if its the right file type
		if fType in fileList[i]:
			fileListNew.append(fileList[i])
	
	#Sorting
	fileListNew.sort()
	fileList=fileListNew
	
	return fileList
				
def lin2winPath(p):
	
	"""Converts Linux Path to Win path (/ -> \\)
	
	Args:
		p (str): Path.
		
	Returns:
		str: Converted path.
		
	"""
	
	r=p.replace("/","\\")
	return r

def win2linPath(p):
	
	"""Converts Win Path to Linux path (\\ -> /)
	
	Args:
		p (str): Path.
		
	Returns:
		str: Converted path.
		
	"""
	
	
	
	r=p.replace("\\","/")
	
	return r

def fixPath(path):
	
	"""Fixes path by expanding user and making sure 
	that path is according to OS definitions.
	
	Args:
		path (str): Path to fix.
		
	Return:
		str: Fixed path.
	
	"""
	
	return os.path.abspath(os.path.expanduser(path))

def leastCommonSubstring(S,T):
	
	"""Find longest common substring.
	
	Taken from http://www.bogotobogo.com/python/python_longest_common_substring_lcs_algorithm_generalized_suffix_tree.php

	Args:
		S (str): String to find substring in.
		T (str): Substring.
		
	Returns:
		str: Least common substring.
		
	"""
	
	m = len(S)
	n = len(T)
	counter = [[0]*(n+1) for x in range(m+1)]
	longest = 0
	lcs_set = set()
	for i in range(m):
		for j in range(n):
			if S[i] == T[j]:
				c = counter[i][j] + 1
				counter[i+1][j+1] = c
				if c > longest:
					lcs_set = set()
					longest = c
					lcs_set.add(S[i-c+1:i+1])
				elif c == longest:
					lcs_set.add(S[i-c+1:i+1])
	lcs_set=list(lcs_set)
	if len(lcs_set)>0:
		lcs_abs=lcs_set[0]
	else:
		lcs_abs=""
		
	return lcs_abs

def str2list(l,dtype="int",openDelim="[",closeDelim="]",sep=","):
	
	"""String with lists/sublists to list filled with integers.
	
	Example:
	
	>>> l="[1,2,[3,4]]"
	>>> str2list(l)
	>>> [1,2,[3,4]]
	
	Args:
		l (str): A list or tuple as a string.
		
	Keyword Args:
		dtype (str): Data type (float,int,str)
		openDelim (str): Opening delimiter.
		closeDelim (str): Closing delimiter.
		sep (str): Seperator between values.
		
	Returns:
		tuple: Tuple containing:
		
			* lnew (list): Converted string.
			* i (int): Last character visited.
		
	"""
	
	#New list
	lnew=[]
	op=False
	j=0
	
	#Find sublists
	for i in range(len(l)):
		k=j+i
		x=l[k]
		
		if x==openDelim:
			
			if op==True:
				#If brackets are already open, recursively start new list
				lrec,jnew=str2list(l[k:],dtype=dtype)	
				
				#Increase j
				j=j+jnew
				
				#Append recursive result
				lnew.append(lrec)
			
			#Start new string
			s=""
			
			#Set open=True
			op=True
			
		elif x==closeDelim:
			#If brackets are closed, append and break
			lnew=appDtype(lnew,s,dtype=dtype)
			op=False
			break
		elif x==sep:
			
			#Append and convert new string
			if len(s)>0:
				lnew=appDtype(lnew,s,dtype=dtype)
				s=""
		else:
			s=s+x
	
	return lnew,i
			
def appDtype(l,s,dtype='int'):
	
	"""Appends string to list and convert to right dtype.
	
	.. note:: Will use ``s.strip()`` before conversion to avoid
	   unnecessary spaces.
	  
	.. note:: Will remove quotes from strings using :py:func:`removeQuoteSignsFromString`.  
	
	Args:
		l (list): A list.
		s (str): String to append.
		
	Keyword Args:
		dtype (str): Data type (float,int,str)
		
	Returns:
		list: List with appended value.
	"""
	
	if s!="":
		s=s.strip()
		if dtype=="int" or dtype==int:
			l.append(int(s))
		elif dtype=="float" or dtype==float:
			l.append(float(s))
		elif dtype=="str" or dtype==str:
			
			# Don't convert strings to strings, that can lead to problems
			s=removeQuoteSignsFromString(s)
			l.append(s)
			
		else:
			printWarning("Not understanding dtype = "+ dtype)
		
	return l	

def removeQuoteSignsFromString(s):
	
	"""Removes Quote signs from string.
	
	Example:
	>>> a='"Test('Test')"'
	>>> removeQuoteSignsFromString()
	>>> 'Test(Test)'
	
	Args:
		s (str): A string.
	
	Returns: 
		str: String without quote signs.
	
	"""
	
	s=s.replace("\"","")
	s=s.replace("\'","")
	
	return s
	
def complValsSimple(l1,l2):
	
	"""Returns complimentary values of two lists.
	
	Args:
		l1 (list): A list.
		l2 (list): Another list.
		
	Returns:
		list: List with complimentary values.
	"""
	
	l=[]
	for i in l1:
		if i not in l2:
			l.append(i)
	
	return l

def complValsFast(l1,l2):
	
	"""Returns complimentary values of two lists, faster version.
	
	Args:
		l1 (list): A list.
		l2 (list): Another list.
		
	Returns:
		list: List with complimentary values.
	"""
	
	matches=matchVals(l1,l2)
	
	for m in matches:
		l1=removeAllOccFromList(l1,m)
	
	return l1

def removeAllOccFromList(l, val):
	
	"""Removes all occurences of value in list.
	
	Args:
		l (list): A list.
		val (value): Value to remove.
		
	Returns:
		list: List without removed values.
	"""
	
	return [value for value in l if value != val]

def matchVals(l1,l2):
	
	"""Returns matching values of two lists.
	
	Args:
		l1 (list): A list.
		l2 (list): Another list.
		
	Returns:
		list: List of matching values.
	"""

	l1=list(l1)
	l2=list(l2)
	
	return list(set(l1).intersection(l2))

def vars2dict(var,loc,filt=[]):	

	"""Builds dict of list of variables (only works if vars are in locals()).
	
	Args:
		var (list): List of variable names.
		loc (dict): Handle to locals().
		
	Returns:
		dict: Built dictionary.
	"""
	
	dic={}
	for name in var:
		dic[name]=loc[name]
		
	return {name: loc[name] for name, val in var.items()}

def objAttr2Dict(obj,attr=[]):
	
	"""Writes all object attributes with names defined in list in the 
	form attributeName = attributeValue into dictionary.
	
	If ``attr=[]``, all attributes are written into dictionary, otherwise
	only the ones specified in ``attr``.
	
	Args:
		obj (object): Object to be printed.
	
	Keyword Args:
		maxL (int): Maximum length threshold.
	
	"""
	
	dic={}
	for item in vars(obj):
		if attr==[]:
			dic[item]=vars(obj)[str(item)]
		else:
			if item in attr:
				dic[item]=vars(obj)[str(item)]
		
		
	return dic


def dict2string(dic,sep="=",newline=False):
	
	"""Build string with variable name and its value from dict.
	
	Args:
		dic (dict): Dictionary.
		
	Keyword Args:
		sep (list): Seperator between variable and value.
		newline (bool): Start newline after each variable.
	
	Returns:
		str: Built string.
	"""
	
	s=""
	
	for i,v in enumerate(dic):
		s=s+v+" "+sep+" "+ str(dic[v])
		if newline:
			if i<len(dic)-1:
				s=s+"\n"
			
	return s		
		
def findIntString(s,idxvec=[],debug=False):			
	
	"""Finds integers in string.
	
	Args:
		s (str): String
		
	Keyword Args:
		idxvec (list): List of already found indices.
		debug (bool): Show debugging output.
	
	Returns:
		list: List of indices if integers in string.
	"""
	
	if idxvec==[]:
		idxvec=list(range(len(s)))
	
	idxs=[]
	found=False
	idx=[]
	for i in idxvec:
		try:
			int(s[i])
			idx.append(i)	
			found=True
		except ValueError:
			if found:
				idxs.append(idx)
				
			found=False
			idx=[]
		
	if found:
		idxs.append(idx)
		
		
	return idxs	

def findDateString(s,sep='',lendate=8,yearreq='',monthreq='',debug=False):
	
	"""Finds date in string. 
	
	Useful for example to find date a filename.
	
	Args:
		s (str): String
		
	Keyword Args:
		sep (str): Separator between date parts.
		lendate (int): Length of dates in characters.
		yearreq (str): Date must contain year, for example ("16").
		monthreq (str): Date must cotain month, for example ("03").
	
	Returns:
		str: Date found, or empty string if date was not found.
		
	"""
	
	idxs=findIntString(s)
	
	for i,idx in enumerate(idxs):
		d=""
		f=False
		if len(sep)>0:
			if i+2<=len(idxs)-1:
				if lenRange(idx)+lenRange(idxs[i+1])+lenRange(idxs[i+2])+3==lendate:
					if s[max(idx):min(idxs[i+1])]==sep and s[max(idx):min(idxs[i+1])]==sep:
						lmin,lmax=rangeLists([idx,idxs[i+1],idxs[i+2]])
						f=True
		else:
			if lenRange(idx)+1==lendate:
				lmin,lmax=rangeLists([idx])
				f=True
		if f:		
			if lmax+1==len(s):
				d=s[lmin:]
			else:
				d=s[lmin:lmax+1]
				
			if yearreq in d and monthreq in d:
				
				return d		
	return ""

	
def lenRange(l):
	
	"""Return the range of values in list
	
	Args:
		l (list): List
	
	Returns:
		float: Range of values of list
		
	"""
	
	return abs(max(l)-min(l))

def rangeLists(ls):
	
	"""Return the range of values in list of lists
	
	Args:
		ls (list): List of lists
	
	Returns:
		tuple: Tuple containing:
		
			* minL (float): Minimum value over all lists
			* maxL (float): Maximum value over all lists
			
	"""
	
	lnew=[]
	for l in ls:
		lnew=lnew+l
	
	minL=min(lnew)
	maxL=max(lnew)
	
	return minL, maxL

def findFn(fn,base,lvlsUp=3,folder=False,debug=False):
	
	"""Finds filename within folder structure
	
	Args:
		fn (str): File name to look for
		base (str): Base path to look in
	
	Keyword Args:
		lvlsUp (int): How many levels to go up
		debug (bool): Debugging flag
		folder (bool): Look for folder
	
	Returns:
		str: Path to preimage
	
	Raises:
		OSError: If file cannot be found
	"""
	
	cwd=os.getcwd()
	
	os.chdir(base)
		
	for i in range(lvlsUp):
		os.chdir('../')
	
	fnFound=""

	for f in os.walk(os.getcwd()):
		f=list(f)
		for k in f:
			if fn in k:
				f.remove([])
				fnFound="".join(f[:-1])+"/"+fn
				if folder:
					if os.path.isdir(fnFound):
						os.chdir(cwd)
						return fnFound
				else:	
					os.chdir(cwd)
					return fnFound
	
	os.chdir(cwd)
	
	if debug:
		print("======= find_folder debugging output =======")
		print("Could not find file ", fn ," . Going to return False")
		print("Closest path found:", fnFound)
		
	return OSError
	
def findPreimage(key,base,lvlsUp=1,fType='tif',debug=False):
	
	"""Finds preimage automatically
	
	Args:
		key (str): Key pattern to look for, e.g. "_pre"
		base (str): Base path to look in
	
	Keyword Args:
		lvlsUp (int): How many levels to go up.
		fType (str): Filetype of preimage
		debug (bool): Debugging flag
	
	Returns:
		str: Path to preimage
	
	Raises:
		OSError: If preimage cannot be found
	"""
	
	folderPre=find_fn('pre',base,lvlsUp=lvlsUp,folder=True,debug=debug)

	if debug:
		print("======= find_preimage debugging output =======")
	
	if not folderPre:
		if debug:
			print("Could not find prefolder with key = ", key ," . Going to return False")
		return folderPre
	else:
		files=get_sorted_folder_list(folderPre,fType)
		if len(files)==0:
			if debug:
				print("Prefolder with key = ", key ," seems to be empty. Going to return False")
			raise OSError("Cannot find preimage")
		else:
			print("Found preimage = ", folderPre+"/"+files[0],  " .")
			return folderPre+"/"+files[0]
		
def updateObj(objBlank,obj,debug=False):
	
	"""Updates object with respect to blank object.
	
	If object does not have attribute, will add 
	attribute from objBlank with value of objBlank.
	
	Args:
		objBlank (object): Object Template
		obj (object): Object to be updated
	
	Keyword Args:
		debug (bool): Print debugging output.
	
	Returns:
		object: Updated object.

	"""
	
	if debug:
		print("======update_obj: updated properties======")
	
	#Going through all attributes blank object
	for item in vars(objBlank):
		
		if not hasattr(obj,str(item)):
			setattr(obj, str(item), vars(objBlank)[str(item)])
			
			if debug:
				print(item, " = ", vars(self)[str(item)])
	
	return obj	
	
def remRepeatsList(l):
	
	"""Removes repeated entries from list. 
	
	Similar to numpy.unique.
	
	Args:
		l (list): List
		
	Returns:
		list: Filtered list.

	"""
	
	return list(set(l))
	 		
def unzipLists(l):
	
	"""Unzips two zipped lists into seperate lists.
	
	Args:
		l (list): Zipped lists
		
	Returns:
		tuple: Tuple containing:
		
		* l1 (list): Unzipped list 1.
		* l2 (list): Unzipped list 2.
		
	"""
	
	l1,l2=list(zip(*l))
	return list(l1),list(l2)
		
def assignIfVal(var,val,valCheck):
	
	"""Assigns val to var if var==valCheck.

	Args:
		var (var): Variable 
		val(any): Value to be assigned
		valCheck(any): Value to be checked for
		
	"""
	
	if var==valCheck:
		return val
	else:
		return var

def enumeratedName(baseName,listOfNames,sep="_"):
	
	"""Generates a new name given a list of names.
	
	Example:
	
	>>> baseName=embryo
	>>> listOfNames=[embryo_1,embryo_5]
	>>> enumeratedName(baseNamem,listOfNames)
	"embryo_6"
	
	Args:
		baseName (str): basename used for naming
		listOfNames(list): List of names
	
	Keyword Args:
		sep (str): Seperator between basename and number
	
	Returns:
		str: New name that has not been used yet
	"""
	
	numbers=[]
	for name in listOfNames:
		if baseName in name and sep in name:
			splittedName=name.strip(baseName).split(sep)
			try:
				numbers.append(int(splittedName[-1]))
			except ValueError:
				pass
	
	if len(numbers)>0:
		newNumber=str(max(numbers)+1)
	else:
		newNumber=str(0)
	
	return baseName+sep+newNumber

def objAttrToList(listOfObjects,AttributeName):
	
	"""Extracts a single object attribute from a list of
	objects and saves it into list.
	
	Args:
		listOfObjects (list): List of objects that all possess the same attribute
		AttributeName (str): Name of attribute to be appended
		
	Returns:
		list: List containing value of attribute of all objects
	"""
	
	l=[]
	for obj in listOfObjects:
		l.append(vars(obj)[str(AttributeName)])
		
	return l

def getAllObjWithAttrVal(listOfObjects,AttributeName,AttributeValue):
	
	"""Filters all objects from a list that have a given attribute value.
	
	Args:
		listOfObjects (list): List of objects that all possess the same attribute.
		AttributeName (str): Name of attribute.
		AttributeValue (str): Value of attribute.
	
	Returns:
		list: List of objects that fulfill requirement.
	"""
	
	vals=objAttrToList(listOfObjects,AttributeName)
	if isinstance(AttributeValue,(list,np.ndarray)):
		b=[]
		for v in vals:
			b.append(compareVectors(v,AttributeValue))
		idx=np.where(b)[0]
	else:			
		idx=np.where(vals==AttributeValue)[0]
	
	return list(np.array(listOfObjects)[idx])

def compareVectors(x,y):
	
	"""Compares two vectors.
	
	Args:
		x (numpy.ndarray): Vector 1.
		y (numpy.ndarray): Vector 2.
		
	Returns:
		bool: True if vectors are identical
		
	"""
	
	return (np.array(x)==np.array(y)).sum()==len(x.flatten())
	
def slashToFn(fn):
	
	"""Append / or \\ to filepath if necessary.
	
	Args:
		fn (str): Filepath 
		
	Returns:
		str: Filepath
	"""
	
	if platform.system() in ["Windows"]:
		s="\\"
	else:
		s="/"
		
	if fn[-1]!=s:
		fn=fn+s
	return fn

def popRange(l,idxStart,idxEnd):
	
	"""Basically ``list.pop()`` for range of indices.
	
	Args:
		l (list): A list
		idxStart (int): Start index of range.
		idxEnd (int): End index of range.
	
	Returns:
		tuple: Tuple containing:
		
			* popped (list): Popped items.
			* l (list): Resulting list.
	
	"""
	
	l=list(l)
	
	if idxStart>idxEnd:
		popped=l[0:idxEnd]
		l[0:idxEnd]=[]
		
		#print popped,l
		
		popped=popped+l[idxStart-idxEnd:]
		l[idxStart-idxEnd:]=[]
	else:
		popped=l[idxStart:idxEnd]
		l[idxStart:idxEnd]=[]
		
	return popped,l
	
def sortListsWithKey(l,keyList):
	
	"""Sorts two lists according to key list.
	
	Example:
	
	>>> l=[1,3,5]
	>>> keyList=[2,5,1]
	
	would result in 
	
	>>> [5,1,3],[1,2,5]
	
	Args:
		l (list): list to be sorted.
		keyList (list): list after which is being sorted.
		
	Returns:
		tuple: Tuple containing:
		
			* sortedList (list): Sorted list
			* sortedKeys (list): Sorted keys
			
	"""
	
	s=sorted(zip(keyList, l), key=lambda keyList:keyList[0])
	
	sortedKeys=[]
	sortedList=[]
	
	for i in range(len(s)):
		sortedKeys.append(s[i][0])
		sortedList.append(s[i][1])
		
	return sortedList,sortedKeys
		
def compareObjAttr(obj1,obj2,):
	
	"""Compare the values of two objects.
	
	Args:
		obj1 (object): First object.
		object2 (object): Second object.
		
	Returns:
		tuple: Tuple containing:
		
			* same (dict): Dictionary of attributes with same values
			* different (dict): Dictionary of attributes with different values
			* notInBoth (dict): Dictionary of attributes that are not in both objects

	"""
	
	same={}
	different={}
	notInBoth={}
	
	for item in vars(obj1):
		
		if item in vars(obj2):
			
			#print item
			
			#Need to check if something is array, then need to use different comparison.
			if isinstance(vars(obj1)[str(item)],(np.ndarray)):
				b=compareArrays(vars(obj1)[str(item)],vars(obj2)[str(item)])	
			else:
				b=(vars(obj1)[str(item)]==vars(obj2)[str(item)])
			
			if b:
				same[str(item)]=[vars(obj1)[str(item)],vars(obj2)[str(item)]]
			else:
				different[str(item)]=[vars(obj1)[str(item)],vars(obj2)[str(item)]]
		else:
			notInBoth[str(item)]=1
	
	for item in vars(obj2):
		if item not in vars(obj1):
			notInBoth[str(item)]=2
			
	return same,different,notInBoth

def compareArrays(arr1,arr2):
	
	"""Converts two lists/arrays into numpy arrays and compares
	then elementwise.
	
	Args:
		arr1 (numpy.ndarray): Some array.
		arr2 (numpy.ndarray): Other arrayy.
		
	Returns:
		bool: True if elemet are the same element wise.
		
	"""
	
	arr1=np.asarray(arr1)
	arr2=np.asarray(arr2)
	
	try:
		return (arr1==arr2).all()
	except:
		(arr1==arr2)
		
def compareROIs(emb1,emb2,byName=True):
	
	"""Compares the list of ROIs between to 
	:py:class:`pyfrp.subclasses.pyfrp_embryo.embryo`` objects.
	
	Args:
		emb1 (pyfrp.subclasses.pyfrp_embryo.embryo): First embryo.
		emb2 (pyfrp.subclasses.pyfrp_embryo.embryo): Second embryo.
		
	Returns:
		tuple: Tuple containing:
		
			* sameAll (list): List of same output of ``compareObjAttr(ROI,ROI2)`` per ROI, see also :py:func:`compareObjAttr`.
			* differentAll (list): List of different output of ``compareObjAttr(ROI,ROI2)`` per ROI, see also :py:func:`compareObjAttr`.
			* notInBothAll (list): List of notInBoth output of ``compareObjAttr(ROI,ROI2)`` per ROI, see also :py:func:`compareObjAttr`.
			* notFound (dict): Dictionary of ROI names that are in ``emb1`` but not ``emb2``.
			
	
	"""
	
	notFound={}
	
	sameAll=[]
	differentAll=[]
	notInBothAll=[]
	
	
	for ROI in emb1.ROIs:
		if byName:
			ROI2=emb2.getROIByName(ROI.name)
		else:
			ROI2=emb2.getROIById(ROI.Id)
		
		if ROI2==None:
			notFound[str(ROI.name)]=1
			continue
		
		same,different,notInBoth=compareObjAttr(ROI,ROI2)
		
		sameAll.append(same)
		differentAll.append(different)
		notInBothAll.append(notInBoth)
		
	return sameAll,differentAll,notInBothAll, notFound

def simpleHist(x,y,bins):

	"""Performs a simple histogram onto x-array.
	
	Args:
		x (numpy.ndarray): x coordinates of data
		y (numpy.ndarray): y coordinates of data
		bins (int):  number of bins
	
	Returns:
		tuple: Tuple containing:
		
			* xBin (numpy.ndarray): Center of bins
			* yBin (numpy.ndarray): Bin Values
		

	"""
	

	xBin=np.linspace(min(x),max(x),bins+1)
			
	iLast=0
	j=1
	
	yBin=[]
	
	for i,xv in enumerate(x):
		try:
			if xv>=xBin[j]:
				yBin.append(np.mean(y[iLast:i]))
				iLast=i
				j=j+1
		except IndexError:
			pass
		
	xBin=np.diff(xBin)+xBin[:-1]
	
	return xBin,np.asarray(yBin)
	
def getMeshfilesDir():
	return 	getModulesDir().replace("modules","meshfiles")

def getSubclassesDir():
	return getModulesDir().replace("modules","subclasses")
	
def getGUIDir():
	return getModulesDir().replace("modules","gui")
	
def getConfDir():
	return getModulesDir().replace("modules","configurations")

def getModulesDir():
	modulePath=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
	return slashToFn(modulePath)

def getMacroDir():
	return slashToFn(getConfDir()+'macros')

def getPathFile():
	return getConfDir()+'paths'

def getFijiBin(fnPath=None):
	return getPath('fijiBin',fnPath=fnPath)

def getOpenscadBin(fnPath=None):
	return getPath('openscadBin',fnPath=fnPath)

def checkIfGmshBin(fn):
	return not os.system(fn + '-1')

def addPathToWinPATHs(path):
	
	"""Adds a path to Windows' PATH list. 
	
	.. note:: Only adds path if file exits.
	
	.. note:: You will need to restart the terminal to 
	   be sure that the change has any effect.
	
	Args:
		path (str): Path to be added.
		
	Returns:
		bool: True if successful.
	"""
	
	if platform.system() not in ["Windows"]:
		printError("OS is not Windows")
		return False
	
	if path in os.environ['PATH']:
		printError("Path is already in PATH.")
		return False
	
	if os.path.exists(path):
		os.system("set PATH=%PATH%;"+path)
		return True
	else:
		printError(path + " does not exist.")
		return False

def addToPATHs(fn):
	
	"""Adds a path to PATH list. 
	
	.. note:: Only adds path if file/folder exits.
	
	Args:
		fn (str): Path to be added.
		
	Returns:
		bool: True if successful.
	"""
	
	# Make sure it does not end on a .
	if fn.endswith('.'):
		fn=fn[:-1]
	
	# Check if folder actually exists
	if not os.path.isdir(fn):
		printError("gmshBin as defined in path file doesn't exist.")
		return False
	
	# Check if it is already in path
	if fn in os.environ['PATH']:
		printNote("gmshBin already in PATH, not going to add it another time.")
		return True
	
	# Append to path
	if platform.system() not in ["Windows"]:
		os.environ['PATH']=os.environ['PATH']+":"+fn
	else:
		os.environ['PATH']=os.environ['PATH']+";"+fn
	
	return True
	
def getPath(identifier,fnPath=None,defaultOutput=""):
	
	"""Extracts path with identifier from path definition file.
	
	If not defined diferently, will first look in configurations/paths,
	then configurations/paths.default.
	
	Args:
		identifier (str): Identifier of path
		
	Keyword Args:
		fnPath (str): Path to path definition file
			
	Returns:
		str: Path

	"""
	
	if fnPath==None:
		fnPath=getPathFile()
	else:
		if not os.path.isfile(fnPath):
			printWarning(fnPath+" does not exist. Will continue with paths defined in default paths files.")
			fnPath=getPathFile()
		
	path=None
	
	with  open (fnPath,'r') as f:
		for line in f:
			if line.strip().startswith(identifier):
				ident,path=line.split('=')
				path=path.strip()
				break
		
	if path==None:
		printWarning("There is no line starting with "+ identifier+"= in "+ fnPath+ ".")
		fnPath=getPathFile()+'.default'
		path=getPath(identifier,fnPath=fnPath)
		
	path=os.path.expanduser(path)
	
	return path

def setPath(identifier,val,fnPath=None):
	
	"""Sets path in path file.
	
	If ``fnPath`` is not given, will use the return of ``getPathFile``.
	
	Args:
		identifier (str): Identifier of path.
		val (str): Value of path.
		
	Keyword Args:
		fnPath (str): Path to path file.
		
	"""
	
	if fnPath==None:
		fnPath=getPathFile()
	else:
		if not os.path.isfile(fnPath):
			printWarning(fnPath+" does not exist. Will continue with paths defined in default paths files.")
			fnPath=getPathFile()
		
	txtLineReplace(fnPath,identifier,identifier+"="+str(val)+"\n")

def printPaths(fnPath=None):
	
	"""Prints out path file.
	
	If ``fnPath`` is not given, will use the return of ``getPathFile``.
		
	Keyword Args:
		fnPath (str): Path to path file.
		
	"""
	
	if fnPath==None:
		fnPath=getPathFile()
		
	with open(fnPath,'rb') as f:
		
		for line in f:
			print(line)
	
def checkPaths(fnPath=None):

	"""Checks if all paths in paths file exist.
	
	If ``fnPath`` is not given, will use the return of ``getPathFile``.
		
	Keyword Args:
		fnPath (str): Path to path file.
		
	"""
	
	if fnPath==None:
		fnPath=getPathFile()
		
	with open(fnPath,'rb') as f:
		
		for line in f:
			try:
				ident,val=line.split("=")
				print(ident, val.strip(), os.path.isfile(val.strip()))
			except ValueError:
				pass
			
def txtLineReplace(filePath, pattern, subst):
		
	"""Replaces line in file that starts with ``pattern`` and substitutes it 
	with ``subst``.
	
	.. note:: Will create temporary file using ``tempfile.mkstemp()``. You should have 
	   read/write access to whereever ``mkstemp`` is putting files.
	
	Args:
		filePath (str): Filename.
		pattern (str): Pattern to be looked for.
		subst (str): String used as a replacement.
			
	"""
	
	
	#Create temp file
	fh, absPath = mkstemp()
	newFile = open(absPath,'w')
	oldFile = open(filePath)
	
	#Loop through file and replace line 
	for line in oldFile:
		
		if line.startswith(pattern):
			newFile.write(line.replace(line, subst))
		else:
			newFile.write(line)
			
	#close temp file
	newFile.close()
	os.close(fh)
	oldFile.close()
		
	#Remove original file
	os.remove(filePath)
	
	#Move new file
	shutil.move(absPath, filePath)
	return
	
def buildEmbryoWizard(fn,ftype,name,nChannel=0,enc="uint16",fnDest=None,createEmbryo=True,recoverIdent=['recover','post'],bleachIdent=['bleach'],preIdent=['pre'],cleanUp=True,colorPrefix='_c00'):
	
	"""Creates embryo object ready for analysis from microscope data.
	
	(1) Extracts microscope data into .tif files
	(2) Builds folder structure
	(3) Moves image files in proper folders
	(4) Creates embryo object and automatically sets filepaths properly
	
	Args:
		fn (str): Path to embryo folder
		ftype (str): Type of microscopy file, such as lsm or czi
		name (str): Name of embryo
	
	Keyword Args:
		nChannel (int): Defines which channel of the images contains relevant data
		fnDest (str): Path of embryo data structure
		createEmbryo (boo): Flag if embryo object should be created
		recoverIdent (list): List of identifiers for recovery data
		bleachIdent (list): List of identifiers for bleach data
		preIdent (list): List of identifiers for pre-bleach data
		cleanUp (bool): Clean up .tif files from other channels afterwards.
	
		
	Returns:
		pyfrp.subclasses.pyfrp_embryo.embryo: Created Embryo in case of success, otherwise -1

	"""
	
	# Fix path
	fn=slashToFn(fixPath(fn))
	
	# Make sure that input directory exists
	if not os.path.isdir(fn):
		printError(fn+ " does not exist.")
		return -1
	
	# Get all files of right type
	l=glob.glob(fn+'/*'+ftype)
	if len(l)==0:
		printError(fn+ "does not contain images of type " + ftype)
		return -1
	
	# Extract images using Bioformats
	for f in l:
		r=pyfrp_img_module.extractBioFormats(f,fn,debug=True,series=0,channel=nChannel,enc=enc,outputformat='tif')
	if r==-1:
		return -1
	
	# Check for output directory
	if fnDest==None:
		fnDest=fn
	fnDest=slashToFn(fixPath(fnDest))
	
	# Create folder structure
	makeEmbryoFolderStruct(fnDest)
	
	
	
	# Sort image files in structure
	sortImageFiles(fn,fnDest,ftype,nChannel=nChannel,recoverIdent=recoverIdent,bleachIdent=bleachIdent,preIdent=preIdent,colorPrefix=colorPrefix,cleanUp=cleanUp)
	
	if createEmbryo:
		from pyfrp.subclasses import pyfrp_embryo
		
		emb=pyfrp_embryo.embryo(name)
		emb.setDataFolder(fnDest+'recover/')
		a=emb.newAnalysis()
		a.setFnPre(fnDest+'pre/')
		
		return emb
	else:
		return 1
	
def sortImageFiles(fn,fnDest,ftype,recoverIdent=['recover','post'],bleachIdent=['bleach'],preIdent=['pre'],nChannel=1,debug=False,colorPrefix='_c00',cleanUp=True):
	
	"""Sorts all image data in fn into the respective folders of embryo project.
	
	Args:
		fn (str): Path to folder containing images
		fnDest (str): Path of embryo project.
		ftype (str): Type of microscopy file, such as lsm or czi

	Keyword Args:
		recoverIdent (list): List of identifiers for recovery data
		bleachIdent (list): List of identifiers for bleach data
		preIdent (list): List of identifiers for pre-bleach data
		nChannel (int): Defines which channel of the images contains relevant data
		debug (bool): Debugging flag
		colorPrefix (str): Defines how to detect if multichannel or not
		cleanUp (bool): Clean up .tif files from other channels afterwards.
		
	Returns:
		int: 0

	"""
	
	#recoverMulti,preMulti,bleachMulti=checkDataMultiChannel(fn,recoverIdent=recoverIdent,bleachIdent=bleachIdent,preIdent=preIdent,colorPrefix=colorPrefix)
	recoverMulti,preMulti,bleachMulti=False,False,False
	
	
	moveImageFiles(fn,'recover','tif',recoverIdent,recoverMulti,fnDest=fnDest,debug=debug,colorPrefix=colorPrefix)
	moveImageFiles(fn,'pre','tif',preIdent,preMulti,fnDest=fnDest,debug=debug,colorPrefix=colorPrefix)
	moveImageFiles(fn,'bleach','tif',bleachIdent,bleachMulti,fnDest=fnDest,debug=debug,colorPrefix=colorPrefix)
	moveImageFiles(fn,'lsm',ftype,[""],False,fnDest=fnDest,debug=debug,colorPrefix=colorPrefix)
	
	if cleanUp:
		cleanUpImageFiles(fn,'tif',ident=recoverIdent,debug=debug,colorPrefix=colorPrefix)
		cleanUpImageFiles(fn,'tif',ident=bleachIdent,debug=debug,colorPrefix=colorPrefix)
		cleanUpImageFiles(fn,'tif',ident=preIdent,debug=debug,colorPrefix=colorPrefix)
		
	return 

def cleanUpImageFiles(fn,ftype,ident=None,debug=False,colorPrefix='_c00',nChannel=None):
	
	"""Removes all image files fullfilling ``*ident*colorPrefix*ftype`` from ``fn``.
	
	If ``nChannel=None``, will remove all files of ftype.
	
	Args:
		fn (str): Path to folder containing images.
		ftype (str): Type of file, for example "tif".
		

	Keyword Args:
		debug (bool): Debugging flag
		colorPrefix (str): String prefix before channel number.
		nChannel (int): Defines which channel to delete
		ident (list): List of identifiers, for example ["recover","post"].
		
	Returns:
		int: Returns 0 if success, -1 if error

	"""
	
	debugFlag=debug*' -v '

	if nChannel!=None:
		colorFlag='*'+colorPrefix+str(nChannel)
	else:
		colorFlag=''
		
	rs=[1]	
	if ident!=None:
		for ind in ident:
			removeListOfFiles(glob.glob(fn+'*'+ind+colorFlag+'*.'+ftype))
					
	else:
		removeListOfFiles(glob.glob(fn+'*'+colorFlag+'*.'+ftype))
		
	return min(rs)
	
	
def moveImageFiles(fn,fnTarget,ftype,ident,isMulti,fnDest=None,debug=False,colorPrefix='_c00',nChannel=1):
	
	"""Moves all image files fullfilling *ident*(isMulti*colorPrefix)* to fnDest+fnTarget.
	
	Args:
		fn (str): Path to folder containing images.\n
		fnTarget (str): Name of folder files should go in, for example "recover".
		ftype (str): Type of file, for example "tif".
		ident (list): List of identifiers, for example ["recover","post"].
		isMulti (bool): Flag if images are multichannel or not.

	Keyword Args:
		fnDest (str): Path containing fnTarget
		debug (bool): Debugging flag
		colorPrefix (str): Defines how to detect if multichannel or not
		nChannel (int): Defines which channel of the images contains relevant data
		
	Returns:
		int: Returns 0 if success, -1 if error

	"""
	
	if fnDest==None:
		fnDest=fn
	
	debugFlag=debug*' -v '
	rs=[1]
	
	for ind in ident:
			
		colorFlag=isMulti*('*'+colorPrefix+str(nChannel))
		
		#cmd='mv '+debugFlag+fn+'*'+ind+colorFlag+'*.'+ftype +' ' + fnDest + fnTarget+'/'
		
		moveListOfFiles(glob.glob(fn+'*'+ind+colorFlag+'*.'+ftype),slashToFn(fixPath(fnDest+fnTarget)))
		
		#try:
			#os.system(cmd)
			#r=1
		#except:
			#printError("Something went wrong executing: " + cmd)
			#r=0
	
	return min(rs)

def copyListOfFiles(l,dest):
	
	"""Copies list of files using shutil.copy.
	
	Args:
		l (list): List of files.
		dest (str): Destination of files.
		
	Returns:
		int: Returns 0 if success, -1 if error
	"""
	
	r=0
	for f in l:
		try:
			shutil.copy(f,dest)
		except:
			printWarning("copyListOfFiles: Was not able to copy: " + f + " -> " + dest)
			r=-1
			
	return r		

def removeListOfFiles(l):
	
	"""Removes list of files using os.remove.
	
	Args:
		l (list): List of files.
		
	Returns:
		int: Returns 0 if success, -1 if error
	"""
	
	r=0
	for f in l:
		try:
			os.remove(f)
		except:
			printWarning("removeListOfFiles: Was not able to remove: " + f )
			r=-1
			
	return r		

def moveListOfFiles(l,dest):
	
	"""Moves list of files using shutil.move.
	
	Args:
		l (list): List of files.
		dest (str): Destination of files.
		
	Returns:
		int: Returns 0 if success, -1 if error
	"""
	
	r=0
	for f in l:
		try:
			shutil.move(f,dest)
		except:
			printWarning("moveListOfFiles: Was not able to move: " + f + " -> " + dest)
			r=-1
			
	return r		
			
def checkDataMultiChannel(fn,recoverIdent=['recover','post'],bleachIdent=['bleach'],preIdent=['pre'],colorPrefix='_c00'):
	
	"""Checks if extracted bleach/pre/recover microscopy data are multichannel.
	
	Args:
		fn (str): Path to folder containing images

	Keyword Args:
		recoverIdent (list): List of identifiers for recovery data
		bleachIdent (list): List of identifiers for bleach data
		preIdent (list): List of identifiers for pre-bleach data
		colorPrefix (str): Defines how to detect if multichannel or not
				
	Returns:
		tuple: Tuple containing:
		
			* recoverMulti (bool): True if recover is multichannel
			* preMulti (bool): True if pre is multichannel
			* bleachMulti (bool): True if bleach is multichannel
		
	"""
	
	l=getSortedFileList(fn,'.tif')
	
	recoverMulti=checkMultiChannel(fn,recoverIdent,colorPrefix=colorPrefix)
	preMulti=checkMultiChannel(fn,preIdent,colorPrefix=colorPrefix)
	bleachMulti=checkMultiChannel(fn,bleachIdent,colorPrefix=colorPrefix)
	
	return recoverMulti,preMulti,bleachMulti
	
	

def checkMultiChannel(fn,ident,colorPrefix='_c00'):
	
	"""Checks if extracted microscopy with identifier are multichannel.
	
	Args:
		fn (str): Path to folder containing images
		ident (list): List of identifiers, for example ["recover","post"]
	
	Keyword Args:
		colorPrefix (str): Defines how to detect if multichannel or not
		
	Returns:
		bool: True if multichannel, False if not

	"""
	
	l=getSortedFileList(fn,'.tif')
	for f in l:
		for i in ident:
			if colorPrefix in f:
				return True
	return False		
			

def makeEmbryoFolderStruct(fn):
	
	"""Creates default folder structure for embryo object.
	
	 fn \n
	 \|--recover \n
         \|--pre \n
	 \|--bleach \n
	 \|--lsm \n
	 \|--meshfiles \n
	
	Args:
		fn (str): Path to embryo folder
			
	Returns:
		int: 0

	"""
	
	mkdir(fn)
	mkdir(fn+"recover")
	mkdir(fn+"pre")
	mkdir(fn+"bleach")
	mkdir(fn+"meshfiles")
	mkdir(fn+"lsm")
	
	return 0
	
def translateNPFloat(x):
	
	"""Translates string into numpy float.
	
	If string==+/-'inf', will return +/- numpy.inf,
	otherwise will leave x unchanged. \n
	This function is necessary to prevent Sphinx from
	not being able to import modules because np.inf is an 
	default value for an keyword arg, but numpy is mocked.
	
	Args:
		x (float): Input number/string.
			
	Returns:
		float: Translated float.

	"""
	
	return np.float(x)
		
def getNLargest(x,N):
	
	"""Returns N largest values in array/list.
	
	Args:
		x (numpy.nparray): An array.
		N (int): Number of values.
		
	Returns:
		list: List containing N largest numbers.
	"""
	
	x=list(x)
	
	return sorted(x)[-N:]

def getIdxOfNLargest(x,N):
	
	"""Returns indices of N largest values in array/list.
	
	Args:
		x (numpy.nparray): An array.
		N (int): Number of values.
		
	Returns:
		tuple: Tuple containing:
		
			* list: List containing N largest numbers.
			* list: List containing indices of N largest numbers.
			
	"""
	
	x=list(x)
	nLargest=getNLargest(x,N)
	
	indices=[]
	for n in nLargest:
		indices.append(x.index(n))
		
		#Need to fill current index with None, in case there are
		#some value is multiple times in x
		x[indices[-1]]=None
		
	return nLargest, indices

def modIdx(i,l):
	
	"""Returns index of list when 
	input is larger than list by returning the modulo of the length of
	the list. 
	
	Useful if lists refer to loop etc.
	
	Args:
		i (int): Index.
		l (list): Some list.
		
	Returns:
		int: New index.
	
	"""
	
	return np.mod(i,len(l))

def mkdir(fn):
	
	"""Tries to make folder if not already existent.
	
	Args:
		fn (str): Path of folder to create.
		
	Returns:
		bool: True if success, False otherwise.
		
	"""
	
	try:
		os.mkdir(fn)
		return True
	except OSError:
		printWarning("Could not create folder "+fn+" .")
		return False

def listToDtype(row,dtype):
	
	"""Converts all entries in list to specified datatype.
	
	Args:
		row (list): Some list.
		dtype (dtype): Some datatype.
		
	Returns:
		list: Converted list.
	"""
	
	rowNew=[]
	for r in row:
		rowNew.append(dtype(r))
	return rowNew		
