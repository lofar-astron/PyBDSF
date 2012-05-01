# $Id$

# Add current directory to -I path
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(_cbdsm MODULE 
  Fitter_dn2g.cc
  Fitter_dnsg.cc
  Fitter_lmder.cc
  MGFunction1.cc
  MGFunction2.cc
  cbdsm_main.cc
  stat.cc
  num_util/num_util.cpp)

target_link_libraries(_cbdsm
  ${Boost_PYTHON_LIBRARY}
  minpack
  port3
)

set_target_properties(_cbdsm PROPERTIES
  PREFIX ""
  LINKER_LANGUAGE Fortran)

if (APPLE)
  set_target_properties(_cbdsm PROPERTIES
    LINK_FLAGS "-undefined dynamic_lookup")
endif (APPLE)   

install(TARGETS _cbdsm
  DESTINATION ${PYTHON_INSTALL_DIR}/lofar/bdsm)


