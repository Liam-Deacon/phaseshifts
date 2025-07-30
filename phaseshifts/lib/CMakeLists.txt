cmake_minimum_required(VERSION 3.15)
project(libphsh LANGUAGES Fortran)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

# Build the Fortran extension as a Python module
Python_add_library(libphsh MODULE libphsh.f WITH_SOABI)
install(TARGETS libphsh DESTINATION phaseshifts/lib)
