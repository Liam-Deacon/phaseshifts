# cmake_minimum_required(VERSION 3.15)
# project(libphsh LANGUAGES Fortran)

# find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

# # Build the Fortran extension as a Python module
# Python_add_library(libphsh MODULE libphsh.f WITH_SOABI)
# install(TARGETS libphsh DESTINATION phaseshifts/lib)

include(FortranCInterface)
FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")

add_library(libphsh MODULE libphsh.f)
set_target_properties(libphsh PROPERTIES
    PREFIX ""
    OUTPUT_NAME "libphsh"
    SUFFIX "${PYTHON_EXTENSION_MODULE_SUFFIX}")
install(
    TARGETS libphsh
    COMPONENT Runtime
    # DESTINATION "${SKBUILD_LIB_DIR}"  # <= previous location
    LIBRARY DESTINATION "phaseshifts/lib"
    RUNTIME DESTINATION "phaseshifts/lib"
)
