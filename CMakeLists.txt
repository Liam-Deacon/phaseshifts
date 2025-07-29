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

# OpenMP option (default OFF due to EXIT statement conflicts in Fortran code)
option(USE_OPENMP "Enable OpenMP support" OFF)

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
set(f2py_module_c "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}module.c")
set(f2py_wrapper_f "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}-f2pywrappers.f")

# Check if source file exists
if(NOT EXISTS "${fortran_src_file}")
  message(FATAL_ERROR "Fortran source file not found: ${fortran_src_file}")
endif()

message(STATUS "Fortran source exists: ${fortran_src_file}")

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
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating f2py wrappers for ${f2py_module_name}"
  VERBATIM
)

# Ensure f2py files are generated before building the module
add_custom_target(generate_f2py_files
  DEPENDS "${f2py_module_c}" "${f2py_wrapper_f}"
  COMMENT "Ensuring f2py files are generated"
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

# Platform-specific linking and compiler flags
if (UNIX)
  if (APPLE)
    set_target_properties("${f2py_module_name}" PROPERTIES
    LINK_FLAGS  '-Wl,-dylib,-undefined,dynamic_lookup')
  else()
    # Linux-specific flags
    if(USE_OPENMP)
      find_package(OpenMP QUIET)
      if(OpenMP_Fortran_FOUND)
        target_link_libraries("${f2py_module_name}" PUBLIC OpenMP::OpenMP_Fortran)
        message(STATUS "OpenMP found and linked")
      else()
        # Fallback: try to link gomp directly if available
        find_library(GOMP_LIBRARY gomp)
        if(GOMP_LIBRARY)
          target_link_libraries("${f2py_module_name}" PUBLIC ${GOMP_LIBRARY})
          target_compile_options("${f2py_module_name}" PRIVATE -fopenmp)
          message(STATUS "Using gomp library: ${GOMP_LIBRARY}")
        else()
          message(WARNING "OpenMP requested but not found - building without OpenMP support")
        endif()
      endif()
    else()
      message(STATUS "OpenMP disabled (USE_OPENMP=OFF)")
    endif()
    set_target_properties("${f2py_module_name}" PROPERTIES
      LINK_FLAGS  '-Wl,--allow-shlib-undefined')
  endif()
endif()

add_dependencies("${f2py_module_name}" genpyf generate_f2py_files)

# Debug output to help diagnose build issues
message(STATUS "Building target: ${f2py_module_name}")
message(STATUS "Fortran source: ${fortran_src_file}")
message(STATUS "F2PY module C: ${f2py_module_c}")
message(STATUS "F2PY wrapper F: ${f2py_wrapper_f}")
message(STATUS "Python extension suffix: ${PYTHON_EXTENSION_MODULE_SUFFIX}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "OpenMP support: ${USE_OPENMP}")
if(NOT USE_OPENMP)
  message(STATUS "Note: To enable OpenMP, use -DUSE_OPENMP=ON (may cause build failures due to EXIT statements in Fortran code)")
endif()

# --- phshift2007 binary build and packaging ---
option(ENABLE_PHSHIFT2007_BINARIES "Build phshift2007 binaries" ON)

if(ENABLE_PHSHIFT2007_BINARIES)
  set(PHSH2007_URL "https://www.icts.hkbu.edu.hk/VanHove_files/leed/phshift2007.zip")
  set(PHSH2007_ZIP "${CMAKE_BINARY_DIR}/phshift2007.zip")
  set(PHSH2007_SRC_DIR "${CMAKE_BINARY_DIR}/phshift2007_src")
  set(PHSH2007_BUILD_DIR "${CMAKE_BINARY_DIR}/phshift2007_build")

  # Download phshift2007.zip if not present
  if(NOT EXISTS "${PHSH2007_ZIP}")
      file(DOWNLOAD "${PHSH2007_URL}" "${PHSH2007_ZIP}")
  endif()

  # Extract phshift2007.zip
  file(MAKE_DIRECTORY "${PHSH2007_SRC_DIR}")
  file(ARCHIVE_EXTRACT INPUT "${PHSH2007_ZIP}" DESTINATION "${PHSH2007_SRC_DIR}")

  # Split .ab3 and .ab4 files into .for files using Python script
  set(SPLIT_SCRIPT "${CMAKE_SOURCE_DIR}/phaseshifts/lib/helpers/split_ab.py")
  add_custom_command(
      OUTPUT "${PHSH2007_BUILD_DIR}/phsh0.for" "${PHSH2007_BUILD_DIR}/phsh1.for" "${PHSH2007_BUILD_DIR}/phsh2cav.for" "${PHSH2007_BUILD_DIR}/phsh2rel.for" "${PHSH2007_BUILD_DIR}/phsh2wil.for" "${PHSH2007_BUILD_DIR}/phsh3.for"
      COMMAND ${Python_EXECUTABLE} "${SPLIT_SCRIPT}" "${PHSH2007_SRC_DIR}/psprog.ab3" "${PHSH2007_BUILD_DIR}"
      COMMAND ${Python_EXECUTABLE} "${SPLIT_SCRIPT}" "${PHSH2007_SRC_DIR}/psprog.ab4" "${PHSH2007_BUILD_DIR}"
      DEPENDS "${PHSH2007_SRC_DIR}/psprog.ab3" "${PHSH2007_SRC_DIR}/psprog.ab4" "${SPLIT_SCRIPT}"
      COMMENT "Splitting .ab3/.ab4 files into .for files for phshift2007"
      VERBATIM
  )
  add_custom_target(phsh_sources ALL DEPENDS 
      "${PHSH2007_BUILD_DIR}/phsh0.for" "${PHSH2007_BUILD_DIR}/phsh1.for" "${PHSH2007_BUILD_DIR}/phsh2cav.for" "${PHSH2007_BUILD_DIR}/phsh2rel.for" "${PHSH2007_BUILD_DIR}/phsh2wil.for" "${PHSH2007_BUILD_DIR}/phsh3.for"
  )

  # --- phshift2007 binary build ---
  set(PHSH_BIN_DIR "${CMAKE_BINARY_DIR}/bin")
  file(MAKE_DIRECTORY "${PHSH_BIN_DIR}")

  set(PHSH_BINARIES phsh0 phsh1 phsh2wil phsh2cav phsh2rel phsh3)
  set(PHSH_SOURCES
      "${PHSH2007_BUILD_DIR}/phsh0.for"
      "${PHSH2007_BUILD_DIR}/phsh1.for"
      "${PHSH2007_BUILD_DIR}/phsh2wil.for"
      "${PHSH2007_BUILD_DIR}/phsh2cav.for"
      "${PHSH2007_BUILD_DIR}/phsh2rel.for"
      "${PHSH2007_BUILD_DIR}/phsh3.for"
  )

  foreach(bin IN LISTS PHSH_BINARIES)
      add_executable(${bin} "${PHSH2007_BUILD_DIR}/${bin}.for")
      set_target_properties(${bin} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PHSH_BIN_DIR}")
      add_dependencies(${bin} phsh_sources)
      # Platform-specific Fortran flags
      if(MSVC)
          # Windows: minimal flags
          target_compile_options(${bin} PRIVATE /O2)
      elseif(APPLE)
          target_compile_options(${bin} PRIVATE -Wall -std=legacy)
      else()
          target_compile_options(${bin} PRIVATE -Wall -Wno-unused-label -Wno-tabs -Wno-unused-variable -Wno-unused-dummy-argument -fcheck=bounds -frecursive -std=legacy)
          target_link_options(${bin} PRIVATE -static-libgcc -static-libgfortran)
      endif()
  endforeach()

  add_custom_target(phshift2007 ALL DEPENDS ${PHSH_BINARIES})

  # Install phshift2007 binaries to bin/
  install(TARGETS ${PHSH_BINARIES}
          RUNTIME DESTINATION bin
          COMPONENT Runtime)
endif()

# Only install the compiled extension - no data files
install(TARGETS "${f2py_module_name}"
        DESTINATION phaseshifts/lib/
        COMPONENT Runtime)
