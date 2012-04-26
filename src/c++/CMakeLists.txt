# $Id$

file(GLOB cbdsm_sources *.cc)

add_library(_cbdsm MODULE
  ${cbdsm_sources}
)

target_link_libraries(_cbdsm
  minpack
  port3
  ${NUM_UTIL_LIBRARIES}
  ${BOOST_PYTHON_LIBRARY}
)

set_target_properties(_cbdsm PROPERTIES
  PREFIX ""
  LINKER_LANGUAGE Fortran
)

if (APPLE)
  set_target_properties(_cbdsm PROPERTIES
    LINK_FLAGS "-undefined dynamic_lookup")
endif (APPLE)   

install(TARGETS _cbdsm
  DESTINATION lib/python/bdsm
)
