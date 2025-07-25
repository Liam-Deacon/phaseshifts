# Modified from https://numpy.org/doc/stable/f2py/buildtools/skbuild.html#cmake-modules-only

### setup project ###
cmake_minimum_required(VERSION 3.15)

# Windows-specific configuration for mixed toolchains
if(WIN32)
  # Force MinGW Makefiles generator on Windows to avoid MSVC/gfortran mixing
  if(NOT CMAKE_GENERATOR STREQUAL "MinGW Makefiles" AND NOT CMAKE_GENERATOR STREQUAL "MSYS Makefiles")
    message(STATUS "Windows detected - forcing MinGW Makefiles generator for compatibility")
    set(CMAKE_GENERATOR "MinGW Makefiles" CACHE STRING "Generator" FORCE)
  endif()

  # Ensure gfortran is found on Windows - check multiple common paths
  if(NOT CMAKE_Fortran_COMPILER)
    find_program(GFORTRAN_EXECUTABLE
      NAMES gfortran gfortran.exe
      PATHS
        "C:/mingw64/bin"
        "C:/tools/mingw64/bin"
        "C:/msys64/mingw64/bin"
        "C:/Program Files/mingw-w64/x86_64-*/mingw64/bin"
        "C:/msys64/ucrt64/bin"
        "C:/msys64/mingw32/bin"
        "C:/ProgramData/chocolatey/lib/mingw/tools/install/mingw64/bin"
        ENV PATH
      DOC "gfortran compiler"
    )
    if(GFORTRAN_EXECUTABLE)
      set(CMAKE_Fortran_COMPILER "${GFORTRAN_EXECUTABLE}" CACHE FILEPATH "Fortran compiler" FORCE)
      message(STATUS "Found gfortran: ${GFORTRAN_EXECUTABLE}")
    else()
      message(WARNING "gfortran not found. Please install MinGW-w64 or ensure gfortran is in PATH")
    endif()
  endif()

  # Also look for gcc if not already set
  if(NOT CMAKE_C_COMPILER)
    find_program(GCC_EXECUTABLE
      NAMES gcc gcc.exe
      PATHS
        "C:/mingw64/bin"
        "C:/tools/mingw64/bin"
        "C:/msys64/mingw64/bin"
        "C:/Program Files/mingw-w64/x86_64-*/mingw64/bin"
        "C:/msys64/ucrt64/bin"
        "C:/msys64/mingw32/bin"
        "C:/ProgramData/chocolatey/lib/mingw/tools/install/mingw64/bin"
        ENV PATH
      DOC "gcc compiler"
    )
    if(GCC_EXECUTABLE)
      set(CMAKE_C_COMPILER "${GCC_EXECUTABLE}" CACHE FILEPATH "C compiler" FORCE)
      message(STATUS "Found gcc: ${GCC_EXECUTABLE}")
    endif()
  endif()

  # On Windows, prefer all-MSVC or all-MinGW toolchain to avoid linking issues
  if(MSVC AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR
      "Mixed MSVC/gfortran toolchain detected on Windows. "
      "This causes linking failures. Please use one of:\n"
      "1. All-MSVC with Intel Fortran\n"
      "2. All-MinGW (gcc + gfortran)\n"
      "3. Use setup.py instead of CMake for Windows builds\n"
      "Suggestion: Set CMAKE_Fortran_COMPILER to use Intel Fortran or configure MinGW"
    )
  endif()

  # Check for Intel Fortran if using MSVC
  if(MSVC)
    message(STATUS "MSVC detected - checking for compatible Fortran compiler")
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
      message(STATUS "Intel Fortran found - compatible with MSVC")
    else()
      message(WARNING "Non-Intel Fortran with MSVC may cause issues")
    endif()
  endif()
endif()

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
if (WIN32)
  # Windows-specific configuration
  message(STATUS "Configuring for Windows build")

  # Set proper extension for Windows
  set_target_properties("${f2py_module_name}" PROPERTIES SUFFIX ".pyd")

  # Handle mixed toolchain issues
  if(MSVC AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    message(WARNING "Mixed MSVC/gfortran toolchain detected - this may cause linking issues")
    message(STATUS "Consider using Intel Fortran or all-MinGW toolchain")
    # Try to set compatible flags
    set_target_properties("${f2py_module_name}" PROPERTIES
      LINK_FLAGS "/NODEFAULTLIB:MSVCRT")
  endif()

  # Set Windows-specific compiler flags for gfortran
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options("${f2py_module_name}" PRIVATE
      $<$<COMPILE_LANGUAGE:Fortran>:-Wall -frecursive -fcheck=bounds -std=legacy>
    )
    # Note: Avoid -pie, -static-libgcc, -static-libgfortran on Windows
    # as they can cause issues with Python extension modules
  endif()

elseif (UNIX)
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

# Only install the compiled extension - no data files
install(TARGETS "${f2py_module_name}"
        DESTINATION phaseshifts/lib/
        COMPONENT Runtime)
