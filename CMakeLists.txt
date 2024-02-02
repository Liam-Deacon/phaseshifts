# Modified from https://numpy.org/doc/stable/f2py/buildtools/skbuild.html#cmake-modules-only

### setup project ###
cmake_minimum_required(VERSION 3.27)

project(phaseshifts
  VERSION 0.1.8
  DESCRIPTION "LIBPHSH module"
  LANGUAGES C Fortran
  )

enable_language(Fortran)

option(NUMPY_F2PY_INCLUDE_DIR "f2py include directory" "")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/phaseshifts/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

set(Python3_FIND_STRATEGY "LOCATION")
set(Python3_FIND_VIRTUALENV FIRST)
find_package (Python3 COMPONENTS Interpreter Development NumPy  REQUIRED)

message(NOTICE "PYTHON3_EXECUTABLE: ${Python3_EXECUTABLE}")
message(NOTICE "PYTHON_INCLUDE_DIRS: ${Python3_INCLUDE_DIRS}")
message(NOTICE "PYTHON_LIBRARIES: ${Python3_LIBRARIES}")
message(NOTICE "PYTHON_NUMPY_INCLUDE_DIRS: ${Python3_NumPy_INCLUDE_DIRS}")

string(TOLOWER "${CMAKE_Fortran_COMPILER_ID}" FORTRAN_COMPILER_ID)
message(NOTICE "Fortran Compiler Type: ${FORTRAN_COMPILER_ID}")

# Safety net
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n"
  )
endif()

execute_process(
  COMMAND "${Python3_EXECUTABLE}" -c "import platform, sys; print('-'.join([platform.python_implementation().lower(), ''.join(map(str, sys.version_info[:2])), platform.processor().lower(), platform.system().lower()]))"
  OUTPUT_VARIABLE CPYTHON_ABI_TAG
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# F2PY headers
if (NOT NUMPY_F2PY_INCLUDE_DIR)
  execute_process(
    COMMAND "${Python3_EXECUTABLE}"
    -c "import os, numpy.f2py; print(os.path.join(os.path.dirname(numpy.f2py.__file__), 'src') if not hasattr(numpy.f2py, 'get_include') else numpy.f2py.get_include())"
    OUTPUT_VARIABLE F2PY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
else()
  set(F2PY_INCLUDE_DIR "${NUMPY_F2PY_INCLUDE_DIR}")
endif()

message(NOTICE "F2PY_INCLUDE_DIR: ${F2PY_INCLUDE_DIR}")

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
  COMMAND ${Python3_EXECUTABLE} -m "numpy.f2py"
      -m "${f2py_module_name}"
      --lower # Important
      "${fortran_src_file}"
  DEPENDS "${fortran_src_file}" # Fortran source
)

set(PHSHIFT2007_BUILD_COMMAND
${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/phaseshifts/build/phshift2007.py"
)

add_custom_command(
  OUTPUT "bin/phsh0" "bin/phsh1" "bin/phsh2cav" "bin/phsh2rel" "bin/phsh2wil" "bin/phsh3"
  COMMAND ${Python3_EXECUTABLE} -m "phaseshifts.phshift2007"
)

add_custom_target("phshift2007.zip"
  COMMAND ${PHSHIFT2007_BUILD_COMMAND} --action download
  COMMENT "Donwloading phshift2007.zip from Michael Van Hove's LEED Calculations Home Page (if necessary)"
)

add_custom_target("phshift2007"
  COMMAND ${Python3_EXECUTABLE} -m "zipfile" -e "phshift2007.zip" "phshift2007"
  COMMENT "Extracting phshift2007.zip"
  DEPENDS "phshift2007.zip"
)

add_custom_target("phsh"
  COMMAND ${PHSHIFT2007_BUILD_COMMAND} --action extract --output-dir "phshift2007"
  COMMENT "Generating phsh executables"
  BYPRODUCTS
    "phshift2007/phsh0.for"
    "phshift2007/phsh1.for"
    "phshift2007/phsh2cav.for"
    "phshift2007/phsh2rel.for"
    "phshift2007/phsh2wil.for"
    "phshift2007/phsh3.for"
  COMMENT "Reconstructing phsh sources"
  DEPENDS "phshift2007"
)

add_library("${f2py_module_name}" MODULE
            "${f2py_wrapper_f}"
            "${f2py_module_c}"
            "${F2PY_INCLUDE_DIR}/fortranobject.c"
            "${fortran_src_file}")

add_executable("phsh0" "phshift2007/phsh0.for")
add_executable("phsh1" "phshift2007/phsh1.for")
add_executable("phsh2cav" "phshift2007/phsh2cav.for")
add_executable("phsh2rel" "phshift2007/phsh2rel.for")
add_executable("phsh2wil" "phshift2007/phsh2wil.for")
add_executable("phsh3" "phshift2007/phsh3.for")

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC -Wall -Wno-unused-label -Wno-tabs -Wno-unused-variable -Wno-unused-dummy-argument -fcheck=bounds -frecursive -std=legacy")

# set compiler flags
target_compile_options("${f2py_module_name}" PRIVATE
                       $<$<COMPILE_LANGUAGE:Fortran>:-fPIC>
                       $<$<COMPILE_LANGUAGE:C>:-fPIC>)

target_include_directories("${f2py_module_name}" PUBLIC
                           ${F2PY_INCLUDE_DIR}
                           ${Python3_INCLUDE_DIRS}
                           ${Python3_NumPy_INCLUDE_DIRS})

set_target_properties("${f2py_module_name}" PROPERTIES PREFIX "")
# target_link_options("${f2py_module_name}" PRIVATE "-static-libgcc" "-static-libgfortran")

set(LIBPHSH_SUFFIX_PREFIX ".${CPYTHON_ABI_TAG}-${FORTRAN_COMPILER_ID}")

if (UNIX)
  if (APPLE)
    set_target_properties(
      "${f2py_module_name}" PROPERTIES
      LINK_FLAGS  '-Wl,-dylib,-undefined,dynamic_lookup'
      SUFFIX "${LIBPHSH_SUFFIX_PREFIX}.dylib")
  else()  # LINUX
    set_target_properties(
      "${f2py_module_name}" PROPERTIES
      LINK_FLAGS  '-Wl,--allow-shlib-undefined'
      SUFFIX "${LIBPHSH_SUFFIX_PREFIX}.so")
  endif()
else()
  set_target_properties("${f2py_module_name}" PROPERTIES SUFFIX "${LIBPHSH_SUFFIX_PREFIX}.pyd")
endif()

add_dependencies("${f2py_module_name}" genpyf)

# install to source dir
install(TARGETS "${f2py_module_name}" DESTINATION phaseshifts/lib/)
