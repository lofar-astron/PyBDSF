# $Id$

add_f2py_module(_pytesselate
  pytess_simple.f
  pytess_roundness.f
  DESTINATION ${PYTHON_INSTALL_DIR}/lofar/bdsm)

