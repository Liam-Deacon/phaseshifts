@echo off
REM Build script for Windows - compiles libphsh Fortran extension
REM This script is called by cibuildwheel before-build on Windows

echo Building libphsh Fortran extension for Windows...
echo Current directory: %CD%

REM Navigate to lib directory
cd /d "%~dp0phaseshifts\lib"
if errorlevel 1 (
    echo ERROR: Could not change to phaseshifts\lib directory
    exit /b 1
)

echo Working in: %CD%
dir libphsh.f

REM Check if libphsh.f exists
if not exist libphsh.f (
    echo ERROR: libphsh.f not found in %CD%
    exit /b 1
)

REM Check for gfortran
where gfortran >nul 2>&1
if errorlevel 1 (
    echo ERROR: gfortran not found in PATH
    echo Please install MinGW-w64 or ensure gfortran is available
    exit /b 1
)

echo Found gfortran:
gfortran --version

REM Run f2py to build the extension
echo Running f2py...
python -m numpy.f2py libphsh.f -m libphsh -c --f77flags="-frecursive -fcheck=bounds -std=legacy"
if errorlevel 1 (
    echo ERROR: f2py build failed
    exit /b 1
)

echo Build completed successfully!
dir libphsh*

REM Return to original directory
cd /d "%~dp0"
