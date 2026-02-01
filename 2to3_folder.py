import os
import glob
import sys

files=glob.glob(sys.argv[1]+'*.py')

for i,f in enumerate(files):
	
	print("%i/%i: Fixing %s"%(i,len(files),f))

	cmd="2to3 -w %s"%f
	os.system(cmd)


	
