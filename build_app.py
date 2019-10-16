import os
import sys
import platform

print("This system is running", platform.system(), platform.release(), platform.architecture()[0])
arch1=platform.architecture()[0]
arch2=platform.machine()

if platform.system() in ["Darwin"]:
	name="pyfrap_"+"OSX_"+platform.mac_ver()[0]+"_"+arch1
	os.system("pyinstaller -F -y -n "+ name +" -i logo/pyfdap_icon.icns --windowed pyfdap_app.py")
elif platform.system() in ["Linux"]:
	name="pyfrap_"+"Linux_"+platform.release()+"_"+arch1
	#os.system("pyinstaller -F -y -n "+ name +" -i logo/pyfdap_icon.ico --windowed pyfdap_app.py")
	os.system("pyinstaller -F -y -n "+ name +" --windowed "+sys.argv[1])
elif platform.system() in ["Windows"]:
	name="pyfdap_"+"Win_"+platform.release()+"_"+arch1
	if arch1=="64bit":
		os.system("pyinstaller -F -y -i logo/pyfdap_icon.ico -n" + name +" pyfdap_app.py")
	else:
		os.system("pyinstaller -F -y -i logo/pyfdap_icon.ico -n" + name +" pyfdap_app_win32build.py")
	
