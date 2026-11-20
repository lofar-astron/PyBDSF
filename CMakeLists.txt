cmake_minimum_required(VERSION 3.18)
project(PyBDSF C Fortran)

find_package(PythonExtensions REQUIRED)
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

if(SKBUILD)
    message(STATUS "The project is built using scikit-build")
    set(PYTHON_INCLUDE_DIRS
        ${Python_INCLUDE_DIR}
        ${Python_NumPy_INCLUDE_DIRS}
    )
else()
    set(PYTHON_INCLUDE_DIRS ${Python_INCLUDE_DIRS})
endif()

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
