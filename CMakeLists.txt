cmake_minimum_required(VERSION 3.11)
project(PyBDSF C CXX Fortran)
find_package(PythonExtensions REQUIRED)
# include(FindPythonExtensions)
include(UsePythonExtensions)
find_package(Python3 REQUIRED COMPONENTS Development)
# find_package(NumPy REQUIRED)
find_package(F2PY)
find_package(Boost COMPONENTS python numpy)

# include_directories(${Python3_INCLUDE_DIRS})
# include_directories(${Boost_INCLUDE_DIRS})

if(SKBUILD)
  message(STATUS "The project is built using scikit-build")
endif()


# find_package(Python REQUIRED Development NumPy)
# message(STATUS "Python_INCLUDE_DIRS=${Python_INCLUDE_DIRS}")

set(MINPACK_SRC
    src/minpack/lmder.f
    src/minpack/lmpar.f
    src/minpack/qrfac.f
    src/minpack/qrsolv.f
    src/minpack/enorm.f
    src/minpack/dpmpar.f
)
add_library(minpack STATIC ${MINPACK_SRC})
set_target_properties(minpack PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

set(PORT3_SRC
    src/port3/dnsg.f
    src/port3/dn2g.f
    src/port3/drnsg.f
    src/port3/drn2g.f
    src/port3/d1mach.f
    src/port3/da7sst.f
    src/port3/dc7vfn.f
    src/port3/dd7tpr.f
    src/port3/dd7upd.f
    src/port3/df7hes.f
    src/port3/dg7lit.f
    src/port3/dg7qts.f
    src/port3/ditsum.f
    src/port3/divset.f
    src/port3/dl7itv.f
    src/port3/dl7ivm.f
    src/port3/dl7mst.f
    src/port3/dl7nvr.f
    src/port3/dl7sqr.f
    src/port3/dl7srt.f
    src/port3/dl7svn.f
    src/port3/dl7svx.f
    src/port3/dl7tsq.f
    src/port3/dl7tvm.f
    src/port3/dl7vml.f
    src/port3/dn2cvp.f
    src/port3/dn2lrd.f
    src/port3/dn2rdp.f
    src/port3/do7prd.f
    src/port3/dparck.f
    src/port3/dq7apl.f
    src/port3/dq7rad.f
    src/port3/dq7rfh.f
    src/port3/dr7mdc.f
    src/port3/drldst.f
    src/port3/ds7cpr.f
    src/port3/ds7lup.f
    src/port3/ds7lvm.f
    src/port3/dv2axy.f
    src/port3/dv2nrm.f
    src/port3/dv7cpy.f
    src/port3/dv7dfl.f
    src/port3/dv7prm.f
    src/port3/dv7scl.f
    src/port3/dv7scp.f
    src/port3/dv7swp.f
    src/port3/i1mach.f
    src/port3/i7mdcn.f
    src/port3/stopx.f
)
add_library(port3 STATIC ${PORT3_SRC})
set_target_properties(port3 PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

# message(STATUS "NumPy_INCLUDE_DIRS=${NumPy_INCLUDE_DIRS}")
# include_directories(${Numpy_INCLUDE_DIRS})

# set(PYTESSELATE_MODULE_NAME "_pytesselate")
set(PYTESSELATE_SRC
    src/fortran/pytess_simple.f
    src/fortran/pytess_roundness.f
)
set(PYTESSELATE_MODULE_FILE ${CMAKE_CURRENT_BINARY_DIR}/${PYTESSELATE_MODULE_NAME}${PYTHON_EXTENSION_MODULE_SUFFIX})
add_custom_target(${PYTESSELATE_MODULE_NAME} ALL
  DEPENDS ${PYTESSELATE_MODULE_FILE}
)
add_custom_command(
  OUTPUT ${PYTESSELATE_MODULE_FILE}
  COMMAND ${F2PY_EXECUTABLE}
    -m ${PYTESSELATE_MODULE_NAME}
    -c
    ${PYTESSELATE_SRC}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
# add_library(_pytesselate MODULE ${PYTESSELATE_SRC})
# target_link_libraries(_pytesselate ${Boost_LIBRARIES})
# set_target_properties(_pytesselate PROPERTIES PREFIX "")
python_extension_module(_pytesselate)
# add_python_extension(_pytesselate
#     SOURCES ${PYTESSELATE_SRC}
# )
# install(TARGETS _pytesselate LIBRARY DESTINATION bdsf)
install(TARGETS ${PYTESSELATE_MODULE_NAME} LIBRARY DESTINATION bdsf)

set(CBDSM_SRC
    src/c++/Fitter_dn2g.cc
    src/c++/Fitter_dnsg.cc
    src/c++/Fitter_lmder.cc
    src/c++/MGFunction1.cc
    src/c++/MGFunction2.cc
    src/c++/cbdsm_main.cc
    src/c++/stat.cc
    src/c++/num_util/num_util.cpp
)
include_directories(
    ${Python3_INCLUDE_DIRS}
    ${Python3_NumPy_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/c++
)
add_library(_cbdsm MODULE ${CBDSM_SRC})
# target_include_directories(_cbdsm PUBLIC ${Python3_INCLUDE_DIRS})
target_link_libraries(_cbdsm minpack port3 gfortran ${Boost_LIBRARIES})
# set_target_properties(_cbdsm PROPERTIES PREFIX "")
python_extension_module(_cbdsm)
install(TARGETS _cbdsm LIBRARY DESTINATION bdsf)

set(NATGRID_SRC
    natgrid/Src/natgrid.c
    natgrid/Src/natgridd.c
    natgrid/Src/natgridmodule.c
    natgrid/Src/natgrids.c
    natgrid/Src/nncrunch.c
    natgrid/Src/nncrunchd.c
    natgrid/Src/nncrunchs.c
    natgrid/Src/nnerror.c
    natgrid/Src/nnuser.c
    natgrid/Src/nnuserd.c
    natgrid/Src/nnusers.c
)

add_library(natgridmodule MODULE ${NATGRID_SRC})
target_include_directories(natgridmodule PRIVATE natgrid/Include)
python_extension_module(natgridmodule)
install(TARGETS natgridmodule LIBRARY DESTINATION bdsf/nat)
