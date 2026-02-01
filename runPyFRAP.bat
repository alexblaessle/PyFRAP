REM This script sets gmsh environment variable and then starts PyFRAP
REM This only works if PyFRAP was installed with --gmsh option

SET PATH=%PATH%;%cd%\executables\gmsh
python pyfrp/PyFRAP.py
