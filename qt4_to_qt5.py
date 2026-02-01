"""Script to search for PyQt4 specific strings and replace it with PyQt5 version.

"""

import os
import sys
import shutil
from tempfile import mkstemp

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
		
		if pattern in line:
			newFile.write(line.replace(pattern, subst))
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

def getAllImportantFiles(fn,ftype='.py'):
	
	# Get all files
	filesNew=[]
	for root, dirs, files in os.walk(fn):
		for name in files:
			filesNew.append(os.path.join(root,name))
	
	files=filterStrList(filesNew,ftype)
	
	return files
	
def findAllOccurences(fn,s):	
	
	"""Finds all occurences of string s in file."""
	
	lineNo=[]
	fnFound=[]
	lines=[]
	N=0
	
	try:
		with open(fn,'r') as f:
			
			for i,line in enumerate(f):
				if s in line:
					lineNo.append(i)
					fnFound.append(fn)
					lines.append(line)
					N=N+line.count(s)
	except IOError:
		print("Skipped file ", fn)
	return fnFound,lines,lineNo,N
	
def filterStrList(l,s,endswith=True):
	
	"""Finds all entries in list that contain ``s``. """
	
	lnew=[]
	for x in l:
		if endswith:
			if x.endswith(s):
				lnew.append(x)
		else:
			if s in x:
				lnew.append(x)
				
	return lnew		

def replace_signal(lines,signal,val=""):
	
	"""Replaces signal with new one."""
	
	signal_str=", QtCore.SIGNAL('%s(%s)'), "%(signal,val)
	
	lines_new=[]
	for line in lines:
		if signal_str in line:
			
			
			s,line=line.split('.connect(')
			line=line.replace(signal_str,'.%s.connect('%signal)
			
			n_tab=s.count('\t')
			line=n_tab*'\t'+line
			
			lines_new.append(line)
		else:
			lines_new.append(line)
			
	return lines_new
	

def read_file(fn):
	
	lines=[]
	with open(fn,'r') as f:
		for line in f:
			lines.append(line)
			
	return lines		

def write_file(fn,lines):
	
	with open(fn,'w') as f:
		for line in lines:
			f.write(line)


files=getAllImportantFiles(os.getcwd())

for fn in files:
	lines=read_file(fn)
	
	
	lines=replace_signal(lines,'clicked')
	lines=replace_signal(lines,'stateChanged',val='int')
	
	
	write_file(fn,lines)
	

