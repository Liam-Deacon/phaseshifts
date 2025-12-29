# Modified from https://numpy.org/doc/stable/f2py/buildtools/skbuild.html#cmake-modules-only

### setup project ###
cmake_minimum_required(VERSION 3.15)

#
# Windows-specific configuration for finding compilers BEFORE project()
# NOTE: We can only set compiler paths here, not check compiler IDs
#       Compiler ID checks must happen AFTER project() is called
if(WIN32)
  message(STATUS "Windows detected - configuring Fortran compiler search paths")

  # Find gfortran on Windows if not already specified
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

  # Find gcc on Windows if not already specified
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
endif()

project(phaseshifts
  VERSION 0.1.8
  DESCRIPTION "LIBPHSH module"
  LANGUAGES C Fortran
)

# Post-project() Windows toolchain validation
# NOTE: These checks MUST be after project() since CMAKE_Fortran_COMPILER_ID
#       and MSVC are only set after project() enables the languages
if(WIN32)
  message(STATUS "Windows build - CMake generator: ${CMAKE_GENERATOR}")
  message(STATUS "Windows build - C compiler: ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID})")
  message(STATUS "Windows build - Fortran compiler: ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID})")

  # Check for mixed MSVC/gfortran toolchain (problematic)
  if(MSVC AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    message(WARNING
      "Mixed MSVC/gfortran toolchain detected on Windows. "
      "This may cause linking failures. Consider using one of:\n"
      "1. All-MSVC with Intel Fortran\n"
      "2. All-MinGW (gcc + gfortran) with -G 'MinGW Makefiles'\n"
      "3. Use setup.py instead of CMake for Windows builds"
    )
  endif()

  # Check for Intel Fortran if using MSVC
  if(MSVC)
    message(STATUS "MSVC detected - checking for compatible Fortran compiler")
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
      message(STATUS "Intel Fortran found - compatible with MSVC")
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      message(STATUS "GNU Fortran (gfortran) detected with MSVC - may have issues")
    else()
      message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
    endif()
  endif()
endif()

include(FortranCInterface)
FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")

# Safety net
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n"
  )
endif()

# OpenMP option (default OFF due to EXIT statement conflicts in Fortran code)
option(USE_OPENMP "Enable OpenMP support" OFF)

# Fortran coverage (gcov/gcovr)
option(ENABLE_FORTRAN_COVERAGE "Build Fortran sources with GCC coverage instrumentation" OFF)
if(NOT ENABLE_FORTRAN_COVERAGE AND DEFINED ENV{PHASESHIFTS_FORTRAN_COVERAGE})
  if(NOT "$ENV{PHASESHIFTS_FORTRAN_COVERAGE}" STREQUAL "")
    set(ENABLE_FORTRAN_COVERAGE ON)
    message(STATUS "PHASESHIFTS_FORTRAN_COVERAGE detected - enabling Fortran coverage flags")
  endif()
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
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

# Get the proper Python extension module suffix
# When scikit-build-core is driving the build, use SKBUILD_SOABI (more reliable)
# Otherwise, query Python's sysconfig for EXT_SUFFIX
if(DEFINED SKBUILD_SOABI)
  # scikit-build-core provides the correct SOABI
  if(WIN32)
    set(PYTHON_EXTENSION_MODULE_SUFFIX ".${SKBUILD_SOABI}.pyd")
  else()
    set(PYTHON_EXTENSION_MODULE_SUFFIX ".${SKBUILD_SOABI}${CMAKE_SHARED_MODULE_SUFFIX}")
  endif()
  message(STATUS "Using SKBUILD_SOABI: ${SKBUILD_SOABI}")
else()
  # Fallback: query Python directly
  execute_process(
    COMMAND "${Python_EXECUTABLE}"
    -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX') or '.so')"
    OUTPUT_VARIABLE PYTHON_EXTENSION_MODULE_SUFFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
endif()
message(STATUS "Python extension suffix: ${PYTHON_EXTENSION_MODULE_SUFFIX}")

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

# Normalise paths for CMake include directories on Windows
file(TO_CMAKE_PATH "${NumPy_INCLUDE_DIRS}" NumPy_INCLUDE_DIRS)
file(TO_CMAKE_PATH "${F2PY_INCLUDE_DIR}" F2PY_INCLUDE_DIR)

# Prepping the module
set(f2py_module_name "libphsh")
set(f2py_module_dirpath "phaseshifts/lib")
set(fortran_src_file "${CMAKE_SOURCE_DIR}/${f2py_module_dirpath}/${f2py_module_name}.f")
set(f2py_signature_file "${CMAKE_SOURCE_DIR}/${f2py_module_dirpath}/${f2py_module_name}.pyf")
set(f2py_source_args "${fortran_src_file}")
set(f2py_depends "${fortran_src_file}")
if(EXISTS "${f2py_signature_file}")
  list(PREPEND f2py_source_args "${f2py_signature_file}")
  list(APPEND f2py_depends "${f2py_signature_file}")
endif()
set(f2py_generated_c "${CMAKE_SOURCE_DIR}/${f2py_module_dirpath}/${f2py_module_name}module.c")
set(f2py_generated_f "${CMAKE_SOURCE_DIR}/${f2py_module_dirpath}/${f2py_module_name}-f2pywrappers.f")
set(use_pre_generated_f2py FALSE)
if(Python_VERSION VERSION_LESS "3.9" AND EXISTS "${f2py_generated_c}" AND EXISTS "${f2py_generated_f}")
  set(use_pre_generated_f2py TRUE)
  set(f2py_module_c "${f2py_generated_c}")
  set(f2py_wrapper_f "${f2py_generated_f}")
else()
  set(f2py_module_c "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}module.c")
  set(f2py_wrapper_f "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}-f2pywrappers.f")
endif()

# Check if source file exists
if(NOT EXISTS "${fortran_src_file}")
  message(FATAL_ERROR "Fortran source file not found: ${fortran_src_file}")
endif()

message(STATUS "Fortran source exists: ${fortran_src_file}")

# Target for enforcing dependencies
add_custom_target(genpyf
  DEPENDS ${f2py_depends}
)
if(NOT use_pre_generated_f2py)
  add_custom_command(
    OUTPUT "${f2py_module_c}" "${f2py_wrapper_f}"
    COMMAND ${PYTHON_EXECUTABLE} -m "numpy.f2py"
        -m "${f2py_module_name}"
        --lower # Important
        ${f2py_source_args}
    DEPENDS ${f2py_depends} # Fortran/signature source
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating f2py wrappers for ${f2py_module_name}"
    VERBATIM
  )
endif()

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

# --- Windows / MinGW-w64 specifics ------------------------------------------
if(WIN32 AND MINGW)
  # 1.  Make sure MinGW is the linker, not MSVC
  set(CMAKE_SYSTEM_NAME Windows)
  # 2.  Ask gfortran to link its own runtime statically
  target_link_options(${f2py_module_name} PRIVATE
        -static-libgfortran      # Fortran runtime
        -static-libgcc           # GCC support library
        -static                  # avoid external CRT DLLs
  )
  # 3.  Also link quadmath if your code uses REAL*16 or EXP/LN variants
  find_library(QUADMATH_LIB quadmath)
  if(QUADMATH_LIB)
      target_link_libraries(${f2py_module_name} PRIVATE ${QUADMATH_LIB})
  endif()
  # 4.  Export all symbols so Python can find init routine
  set_target_properties(${f2py_module_name} PROPERTIES
        CXX_VISIBILITY_PRESET default
        WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

if(WIN32)
  # On Windows, we need to link against the Python library explicitly
  target_link_libraries("${f2py_module_name}" PRIVATE Python::Module)
endif()

target_include_directories("${f2py_module_name}" PUBLIC
                           ${F2PY_INCLUDE_DIR}
                           ${NumPy_INCLUDE_DIRS}
                           ${Python_INCLUDE_DIRS})

if(WIN32)
  # Windows Python extension modules use .pyd
  set_target_properties("${f2py_module_name}" PROPERTIES SUFFIX ".pyd")
elseif(APPLE)
  # macOS Python extension modules use .so
  set_target_properties("${f2py_module_name}" PROPERTIES SUFFIX ".so")
elseif(UNIX)
  # Linux/Unix Python extension modules use .so
  set_target_properties("${f2py_module_name}" PROPERTIES SUFFIX ".so")
else()
  # Fallback: default to .so
  set_target_properties("${f2py_module_name}" PROPERTIES SUFFIX ".so")
endif()
set_target_properties("${f2py_module_name}" PROPERTIES PREFIX "")



# Platform-specific linking and compiler flags
if (WIN32)
  # Windows-specific configuration
  message(STATUS "Configuring for Windows build")

  # Ensure proper extension for Windows (fallback to .pyd if EXT_SUFFIX was empty)
  if(NOT PYTHON_EXTENSION_MODULE_SUFFIX OR PYTHON_EXTENSION_MODULE_SUFFIX STREQUAL ".so")
    set_target_properties("${f2py_module_name}" PROPERTIES SUFFIX ".pyd")
    message(STATUS "Using fallback Windows extension suffix: .pyd")
  endif()

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
    # Use static linking for MinGW runtime libraries to avoid DLL dependency issues
    # This bundles libgfortran, libgcc, and libquadmath into the extension
    target_link_options("${f2py_module_name}" PRIVATE
      -static-libgfortran
      -static-libgcc
    )
  endif()

  # Link Python library on Windows (required - unlike Unix, Windows needs explicit linking)
  target_link_libraries("${f2py_module_name}" PRIVATE ${Python_LIBRARIES})

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

if(ENABLE_FORTRAN_COVERAGE)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    message(STATUS "Fortran coverage enabled (gfortran --coverage)")
    target_compile_options("${f2py_module_name}" PRIVATE
      $<$<COMPILE_LANGUAGE:Fortran>:-g -O0 -fprofile-arcs -ftest-coverage>
    )
    target_link_options("${f2py_module_name}" PRIVATE --coverage)
  else()
    message(WARNING "Fortran coverage requested but compiler is ${CMAKE_Fortran_COMPILER_ID}; skipping coverage flags")
  endif()
endif()


add_dependencies("${f2py_module_name}" genpyf generate_f2py_files)

# --- Optional on‑the‑fly generation of libphsh (disabled by default) ---
option(GENERATE_LIBPHSH "Generate libphsh Fortran library automatically if missing" OFF)
if(GENERATE_LIBPHSH)
  add_custom_target(generate_libphsh
    COMMAND ${CMAKE_COMMAND} -E echo "TODO: generate libphsh from source"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Generating libphsh Fortran library")
  add_dependencies("${f2py_module_name}" generate_libphsh)
endif()

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
option(ENABLE_PHSHIFT2007_BINARIES "Build phshift2007 binaries" OFF)

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
      if(ENABLE_FORTRAN_COVERAGE AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
          target_compile_options(${bin} PRIVATE -g -O0 -fprofile-arcs -ftest-coverage)
          target_link_options(${bin} PRIVATE --coverage)
      endif()
  endforeach()

  add_custom_target(phshift2007 ALL DEPENDS ${PHSH_BINARIES})

  # Install phshift2007 binaries to bin/
  install(TARGETS ${PHSH_BINARIES}
          RUNTIME DESTINATION bin
          COMPONENT Runtime)
endif()

# Only install the compiled extension - no data files
# wheel.install-dir in pyproject.toml is set to "phaseshifts/lib"
# so DESTINATION . installs into phaseshifts/lib/
# For editable mode (inplace), extension is copied to source tree phaseshifts/lib/
install(TARGETS "${f2py_module_name}"
        DESTINATION .
        COMPONENT Runtime)
