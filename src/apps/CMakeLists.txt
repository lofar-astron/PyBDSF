# $Id$

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/pybdsm.in
  ${CMAKE_CURRENT_BINARY_DIR}/pybdsm @ONLY)

lofar_add_bin_scripts(${CMAKE_CURRENT_BINARY_DIR}/pybdsm)
  
