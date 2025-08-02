# cmake_minimum_required(VERSION 3.15)
# project(libphsh LANGUAGES Fortran)

# find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

# # Build the Fortran extension as a Python module
# Python_add_library(libphsh MODULE libphsh.f WITH_SOABI)
# install(TARGETS libphsh DESTINATION phaseshifts/lib)

add_library(libphsh MODULE libphsh.f)
set_target_properties(libphsh PROPERTIES
    PREFIX ""                         # drop the leading “lib”
    OUTPUT_NAME "libphsh")            # final file: libphsh.cpython-...so
install(
    TARGETS libphsh
    COMPONENT Runtime
    # DESTINATION "${SKBUILD_LIB_DIR}"  # <= previous location
    LIBRARY DESTINATION "${SKBUILD_PLATLIB_DIR}/phaseshifts/lib"
    RUNTIME DESTINATION "${SKBUILD_PLATLIB_DIR}/phaseshifts/lib"
)
