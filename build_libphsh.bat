REM Simple build script for windows command
cd .\phaseshifts\lib
dir
f2py -m libphsh -c
cd ..\..
