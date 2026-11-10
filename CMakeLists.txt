cmake_minimum_required(VERSION 3.18)
project(PyBDSF C Fortran)

find_package(PythonExtensions REQUIRED)
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(F2PY REQUIRED)

if(SKBUILD)
    message(STATUS "The project is built using scikit-build")
    set(PYTHON_INCLUDE_DIRS
        ${Python_INCLUDE_DIR}
        ${Python_NumPy_INCLUDE_DIRS}
    )
else()
    set(PYTHON_INCLUDE_DIRS ${Python_INCLUDE_DIRS})
endif()

## -- F2py module: pytesselate -- ##
set(F2PY_MODULE_NAME "_pytesselate")
set(F2PY_MODULE_FILE ${CMAKE_CURRENT_BINARY_DIR}/${F2PY_MODULE_NAME}${PYTHON_EXTENSION_MODULE_SUFFIX})
set(PYTESSELATE_SRC
    ${CMAKE_SOURCE_DIR}/src/fortran/pytess_roundness.f
    ${CMAKE_SOURCE_DIR}/src/fortran/pytess_simple.f
    ${CMAKE_SOURCE_DIR}/src/fortran/constants.h
    )
add_custom_target(${F2PY_MODULE_NAME} ALL
    DEPENDS ${F2PY_MODULE_FILE}
)
add_custom_command(
    OUTPUT ${F2PY_MODULE_FILE}
    COMMAND ${F2PY_EXECUTABLE} -m ${F2PY_MODULE_NAME} -c ${PYTESSELATE_SRC}
    COMMENT "[F2PY] Building python extension module ${F2PY_MODULE_NAME}"
)
install(FILES ${F2PY_MODULE_FILE} DESTINATION bdsf)

## -- Python module: natgrid -- ##
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
target_compile_definitions(natgridmodule PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
target_include_directories(natgridmodule PRIVATE natgrid/Include)
python_extension_module(natgridmodule)
install(TARGETS natgridmodule LIBRARY DESTINATION bdsf/nat)
