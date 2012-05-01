# $Id$

file(GLOB _natgrid_sources Src/*.c)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Include)

add_library(natgridmodule MODULE ${_natgrid_sources})
set_target_properties(natgridmodule PROPERTIES PREFIX "")

if (APPLE)
  set_target_properties(natgridmodule PROPERTIES
    LINK_FLAGS "-undefined dynamic_lookup")
endif (APPLE)   

install(TARGETS natgridmodule DESTINATION ${PYTHON_INSTALL_DIR}/lofar/bdsm)
install(FILES   Lib/nat.py    DESTINATION ${PYTHON_INSTALL_DIR}/lofar/bdsm)
