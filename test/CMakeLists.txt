# $Id$

include(LofarCTest)

# Create symlinks to some Python modules so that they can be found during
# build/test time.
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${PACKAGE_BINARY_DIR}/src/c++/_cbdsm.so
  ${PYTHON_BUILD_DIR}/lofar/bdsm/_cbdsm.so)

execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${PACKAGE_BINARY_DIR}/src/fortran/_pytesselate.so
  ${PYTHON_BUILD_DIR}/lofar/bdsm/_pytesselate.so)

execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${PACKAGE_BINARY_DIR}/src/natgrid/natgridmodule.so
  ${PYTHON_BUILD_DIR}/lofar/bdsm/natgridmodule.so)

execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${PACKAGE_SOURCE_DIR}/src/natgrid/Lib/nat.py
  ${PYTHON_BUILD_DIR}/lofar/bdsm/nat.py)

lofar_add_test(tbdsm_import DEPENDS _cbdsm _pytesselate natgridmodule)
lofar_add_test(tbdsm_process_image DEPENDS _cbdsm _pytesselate natgridmodule)
