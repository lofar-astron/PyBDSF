
## ------------------------------------------------------------------------------
## Project description

## Name of the project
project (BDSM Fortran C CXX)

## Minimum required version of CMake to configure the project
cmake_minimum_required (VERSION 2.8)

## ==============================================================================
##
##  CMake configuration settings
##
## ==============================================================================

if (NOT LUS_SOURCE_DIR)
  find_path (LUS_SOURCE_DIR devel_common/cmake/LUS_CMakeSettings.cmake
    PATHS 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}/../..
    ${CMAKE_CURRENT_SOURCE_DIR}/../../..
    $ENV{LOFARSOFT}
    )
endif (NOT LUS_SOURCE_DIR)

if (LUS_SOURCE_DIR)
  include (${LUS_SOURCE_DIR}/devel_common/cmake/LUS_CMakeSettings.cmake)
else (LUS_SOURCE_DIR)
  message (FATAL_ERROR "Unable to locate additional CMake scripts!")
endif (LUS_SOURCE_DIR)

## ==============================================================================
##
##  Options
##
## ==============================================================================

option (ANAAMIKA_BUILD_TESTS       "Build the test programs?"               YES )
option (ANAAMIKA_ENABLE_TESTING    "Build the test programs?"               YES )
option (ANAAMIKA_VERBOSE_CONFIGURE "Verbose output during configuration?"   NO  )
option (ANAAMIKA_ENABLE_FBDSM      "Build Fortran-based version of BDSM?"   NO  )
option (ANAAMIKA_ENABLE_PYBDSM     "Build Python-enabled version of BDSM?"  YES )

## ==============================================================================
##
##  External dependencies
##
## ==============================================================================

foreach (anaamikaExternalPackages
    Boost
    CASACORE
    CFITSIO
    NumPy
    NumUtil
    GFortran
    PGPLOT
    Python
    WCSLIB
    X11
    )
  lusFindModules (${anaamikaExternalPackages} ${ANAAMIKA_VERBOSE_CONFIGURE})
endforeach (anaamikaExternalPackages)

## ==============================================================================
##
##  Handle options and external dependencies
##
## ==============================================================================

if (NOT HAVE_BOOST OR NOT BOOST_PYTHON_LIBRARY)
  message (STATUS "[BDSM] Unable to build pybdsm - missing Boost!")
  set (ANAAMIKA_ENABLE_PYBDSM FALSE)
endif (NOT HAVE_BOOST OR NOT BOOST_PYTHON_LIBRARY)

if (NOT HAVE_NUMPY)
  message (STATUS "[BDSM] Unable to build pybdsm - missing NumPy!")
  set (ANAAMIKA_ENABLE_PYBDSM FALSE)
endif (NOT HAVE_NUMPY)

if (NOT HAVE_PYTHON)
  message (STATUS "[BDSM] Unable to build pybdsm - missing Python!")
  set (ANAAMIKA_ENABLE_PYBDSM FALSE)
endif (NOT HAVE_PYTHON)

if (NOT NUM_UTIL_LIBRARIES)
  message (STATUS "[BDSM] Unable to build pybdsm - missing NumUtil!")
  set (ANAAMIKA_ENABLE_PYBDSM FALSE)
endif (NOT NUM_UTIL_LIBRARIES)

if (ANAAMIKA_ENABLE_TESTING)
  include (CTest)
  enable_testing() 
endif (ANAAMIKA_ENABLE_TESTING)

# Set Fortran include directories and compile flags.
include_directories (${CMAKE_SOURCE_DIR}/implement/fBDSM/includes)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC -fno-second-underscore")

## ==============================================================================
##
##  Project subdirectories
##
## ==============================================================================

add_subdirectory (implement)
add_subdirectory (apps)

find_package(Sphinx)
if (HAVE_SPHINX)
  add_subdirectory (doc)
endif (HAVE_SPHINX)

## ==============================================================================
##
##  Configuration summary
##
## ==============================================================================

message (STATUS "------------------------------------------------------------")
message (STATUS "[Anaamika] Configuration summary.")
message (STATUS "------------------------------------------------------------")
message (STATUS " System configuration:"                                      )
message (STATUS " .. Processor type ............. = ${CMAKE_SYSTEM_PROCESSOR}")
message (STATUS " .. CMake executable ........... = ${CMAKE_COMMAND}"         )
message (STATUS " .. CMake version .............. = ${CMAKE_VERSION}"         )
message (STATUS " .. System name ................ = ${CMAKE_SYSTEM}"          )
message (STATUS " .. C++ compiler ............... = ${CMAKE_CXX_COMPILER}"    )
message (STATUS " .. C compiler ................. = ${CMAKE_C_COMPILER}"      )
message (STATUS " .. size(void*) ................ = ${CMAKE_SIZEOF_VOID_P}"   )
message (STATUS " Required external packages:"                                )
message (STATUS " .. Boost++ .................... = ${BOOST_FOUND}"           )
message (STATUS " .. CFITSIO .................... = ${CFITSIO_FOUND}"         )
message (STATUS "    .. CFITSIO include dir       = ${CFITSIO_INCLUDES}"      )
message (STATUS "    .. CFITSIO library           = ${CFITSIO_LIBRARIES}"     )
message (STATUS " .. casacore ................... = ${CASACORE_FOUND}"        )
message (STATUS " .. num_util ................... = ${NUM_UTIL_FOUND}"        )
message (STATUS " .. NumPy ...................... = ${NUMPY_FOUND}"           )
message (STATUS " .. Python ..................... = ${PYTHON_FOUND}"          )
message (STATUS "    .. Python version            = ${PYTHON_VERSION}"        )
message (STATUS "    .. Python interpreter        = ${PYTHON_EXECUTABLE}"     )
message (STATUS "    .. Python library            = ${PYTHON_LIBRARIES}"      )
message (STATUS " .. WCSLIB ..................... = ${WCSLIB_FOUND}"          )
message (STATUS "    .. WCSLIB include dir        = ${WCSLIB_INCLUDES}"       )
message (STATUS "    .. WCSLIB library            = ${WCSLIB_LIBRARIES}"      )
message (STATUS "------------------------------------------------------------")
