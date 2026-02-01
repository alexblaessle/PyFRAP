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

"""Terminal module for PyFRAP toolbox. Provides extra functions for a nicer 
custom output inside a Python/bash terminal.
"""

#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

import colorama
import PyQt5.QtGui as QtGui
import numpy as np
import inspect

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================


def printWarning(txt,showCall=True,idx=2):
	
	"""Prints Warning of the form "WARNING: txt", while warning is rendered yellow.
	
	Args:
		txt (str): Text to be printed:
		
	Keyword Args:	
		showCall (bool): Show function in which print function was called.
		idx (int): Traceback index that is given to :py:func:`pyfrp.modules.pyfrp_term_module.getFunctionCall`.
			
	
	"""

	print((colorama.Fore.YELLOW + "WARNING "+showCall*("("+getFunctionCall(idx)+")")+": ") + colorama.Fore.RESET + txt)

def printError(txt,showCall=True,idx=2):
	
	"""Prints Error of the form "ERROR: txt", while error is rendered red.
	
	Args:
		txt (str): Text to be printed:
		
	Keyword Args:	
		showCall (bool): Show function in which print function was called.
		idx (int): Traceback index that is given to :py:func:`pyfrp.modules.pyfrp_term_module.getFunctionCall`.
			
	
	"""
	
	print((colorama.Fore.RED + "ERROR "+showCall*("("+getFunctionCall(idx)+")")+": ") + colorama.Fore.RESET + txt)

def printNote(txt,showCall=True,idx=2):
	
	"""Prints note of the form "NOTE: txt", while note is rendered green.	
	
	Args:
		txt (str): Text to be printed:
		
	Keyword Args:	
		showCall (bool): Show function in which print function was called.
		idx (int): Traceback index that is given to :py:func:`pyfrp.modules.pyfrp_term_module.getFunctionCall`.
			
	
	"""

	print((colorama.Fore.GREEN + "NOTE "+showCall*("("+getFunctionCall(idx)+")")+": ") + colorama.Fore.RESET + txt)
	
def printDict(dic,maxL=5):
	
	"""Prints all dictionary entries in the form key = value.
	
	If attributes are of type ``list`` or ``numpy.ndarray``, will check if the size
	exceeds threshold. If so, will only print type and dimension of attribute.
	
	Args:
		dic (dict): Dictionary to be printed.
		
	Returns:
		bool: True
	
	"""
	
	for k in list(dic.keys()):
		printAttr(k,dic[k],maxL=maxL)
		
	return True	

def printObjAttr(var,obj):
	
	"""Prints single object attribute in the form attributeName = attributeValue.
	Args:
		var (str): Name of attribute.
		obj (object): Object to be printed.
	
	Returns:
		str: Name of attribute.
	
	"""
	
	
	print(var, " = ", vars(obj)[str(var)])
	return var

def printAllObjAttr(obj,maxL=5):
	
	"""Prints all object attributes in the form attributeName = attributeValue.
	
	If attributes are of type ``list`` or ``numpy.ndarray``, will check if the size
	exceeds threshold. If so, will only print type and dimension of attribute.
	
	Args:
		obj (object): Object to be printed.
	
	Keyword Args:
		maxL (int): Maximum length threshold.
	
	"""
	
	for item in vars(obj):
		printAttr(item,vars(obj)[str(item)],maxL=maxL)
	print()	
	return True

def printAttr(name,attr,maxL=5):

	"""Prints single attribute in the form attributeName = attributeValue.
	
	If attributes are of type ``list`` or ``numpy.ndarray``, will check if the size
	exceeds threshold. If so, will only print type and dimension of attribute.
	
	Args:
		name (str): Name of attribute.
		attr (any): Attribute value.
		
	Keyword Args:
		maxL (int): Maximum length threshold.
	
	"""

	if isinstance(attr,(list)):
		if len(attr)>maxL:
			print(name, " = ", getListDetailsString(attr))
			return True
	elif isinstance(attr,(np.ndarray)):
		if min(attr.shape)>maxL:
			print(name, " = ", getArrayDetailsString(attr))
			return True
		
	print(name, " = ", attr)
		
	return True	


def getListDetailsString(l):
	
	"""Returns string saying "List of length x", where x is the length of the list. 
	
	Args:
		l (list): Some list.
	
	Returns:
		str: Printout of type and length.
	"""
	
	return "List of length " + str(len(l))

def getArrayDetailsString(l):
		
	"""Returns string saying "Array of shape x", where x is the shape of the array. 
	
	Args:
		l (numpy.ndarray): Some array.
	
	Returns:
		str: Printout of type and shape.
	"""
	
	return "Array of shape " + str(l.shape)

def printTable(l,header,col=False):
	
	"""Prints table using tabulate package.
	
	If ``col=True``, columns are given via ``l``, otherwise rows are given.
	
	Args:
		l (list): List of rows or columns.
		header (list): List of headers.
		col (bool): Flag on how rows/columns are given.
	
	Returns:
		tuple: Tuple containing:
			
			* header (list): Header of table.
			* table (list): Table as a list of rows.
	
	"""
	
	if col:
		table=[]
		for i in range(len(l[0])):
			row=[]
			for j in range(len(l)):
				row.append(l[j][i])
			table.append(row)
	else:
		table=l
		
	from tabulate import tabulate
			
	print(tabulate(table,headers=header))
		
	return header, table	

def getFunctionCall(idx=1):
	
	"""Returns the name of function or method that was called idx frames outwards.
	
	.. note:: ``idx=0`` will of course return ``getFunctionCall``.
	
	Args:
		idx (int): Steps outwards.
	
	Returns:
		str: Name of function or method.
	"""
	
	frame = inspect.currentframe()
	try:
		return inspect.getouterframes(frame)[idx][3]
	except IndexError:
		return ""
