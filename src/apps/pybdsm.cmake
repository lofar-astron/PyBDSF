#!/bin/bash
#
# This script simply starts the interactive PyBDSM
# IPython shell.

# Check whether called from CEPI/II. If so, init
# the Python libraries before calling PyBDSM.
if [ $APS_LOCAL ]; then
    . ${APS_LOCAL}/login/alias.bash
    use Pythonlibs
fi

# Add the LUS libraries to the relevant paths.
if [ `uname` == "Darwin" ]; then
    DYLD_FALLBACK_LIBRARY_PATH=@CMAKE_INSTALL_PREFIX@/lib${DYLD_FALLBACK_LIBRARY_PATH:+:${DYLD_FALLBACK_LIBRARY_PATH}}
else
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH:+${LD_LIBRARY_PATH}:}@CMAKE_INSTALL_PREFIX@/lib
fi
PYTHONPATH=${PYTHONPATH:+${PYTHONPATH}:}@CMAKE_INSTALL_PREFIX@/lib/python

# And execute pybdsm.py.
exec @PYTHON_EXECUTABLE@ @CMAKE_INSTALL_PREFIX@/lib/python/bdsm/pybdsm.py
