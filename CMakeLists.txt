# Modified from https://numpy.org/doc/stable/f2py/buildtools/skbuild.html#cmake-modules-only

### setup project ###
cmake_minimum_required(VERSION 3.15)

project(phaseshifts
  VERSION 0.1.8
  DESCRIPTION "LIBPHSH module"
  LANGUAGES C Fortran
  )

# Safety net
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n"
  )
endif()

# Ensure scikit-build modules
if (NOT SKBUILD)
  find_package(PythonInterp 3.12 REQUIRED)
  # If skbuild is not the driver; include its utilities in CMAKE_MODULE_PATH
  execute_process(
    COMMAND "${PYTHON_EXECUTABLE}"
    -c "import os, skbuild; print(os.path.dirname(skbuild.__file__))"
    OUTPUT_VARIABLE SKBLD_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  list(APPEND CMAKE_MODULE_PATH "${SKBLD_DIR}/resources/cmake")
  message(STATUS "Looking in ${SKBLD_DIR}/resources/cmake for CMake modules")
endif()

# Use modern scikit-build-core CMake modules
find_package(Python COMPONENTS Interpreter Development REQUIRED)

# NumPy include directory
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}"
  -c "import numpy; print(numpy.get_include())"
  OUTPUT_VARIABLE NumPy_INCLUDE_DIRS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# F2PY headers
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}"
  -c "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Prepping the module
set(f2py_module_name "libphsh")
set(f2py_module_dirpath "phaseshifts/lib")
set(fortran_src_file "${CMAKE_SOURCE_DIR}/${f2py_module_dirpath}/${f2py_module_name}.f")
set(f2py_module_c "${f2py_module_name}module.c")
set(f2py_wrapper_f "${f2py_module_name}-f2pywrappers.f")

# Target for enforcing dependencies
add_custom_target(genpyf
  DEPENDS "${fortran_src_file}"
)
add_custom_command(
  OUTPUT "${f2py_module_c}" "${f2py_wrapper_f}"
  COMMAND ${PYTHON_EXECUTABLE} -m "numpy.f2py"
      -m "${f2py_module_name}"
      --lower # Important
      "${fortran_src_file}"
  DEPENDS "${fortran_src_file}" # Fortran source
)

add_library("${f2py_module_name}" MODULE
            "${f2py_wrapper_f}"
            "${f2py_module_c}"
            "${F2PY_INCLUDE_DIR}/fortranobject.c"
            "${fortran_src_file}")

target_include_directories("${f2py_module_name}" PUBLIC
                           ${F2PY_INCLUDE_DIR}
                           ${NumPy_INCLUDE_DIRS}
                           ${Python_INCLUDE_DIRS})

set_target_properties("${f2py_module_name}" PROPERTIES SUFFIX "${PYTHON_EXTENSION_MODULE_SUFFIX}")
set_target_properties("${f2py_module_name}" PROPERTIES PREFIX "")



if (UNIX)
  if (APPLE)
    set_target_properties("${f2py_module_name}" PROPERTIES
    LINK_FLAGS  '-Wl,-dylib,-undefined,dynamic_lookup')
  else()
    set_target_properties("${f2py_module_name}" PROPERTIES
  LINK_FLAGS  '-Wl,--allow-shlib-undefined')
  endif()
endif()

add_dependencies("${f2py_module_name}" genpyf)

install(TARGETS "${f2py_module_name}" DESTINATION phaseshifts/lib/)
