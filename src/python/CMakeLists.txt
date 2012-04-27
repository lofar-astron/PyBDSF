# $Id$

include(PythonInstall)

python_install(
  __init__.py
  _version.py
  cleanup.py
  collapse.py
  compute_shapelet_coeff.py
  const.py
  fbdsm_fitstable_stuff.py
  functions.py
  gaul2srl.py
  gausfit.py
  image.py
  interface.py
  islands.py
  libs.py
  make_residimage.py
  misc.py
  mylogger.py
  opts.py
  output.py
  output_fbdsm_files.py
  plotresults.py
  polarisation.py
  preprocess.py
  psf_vary.py
  pybdsm.py
  readimage.py
  rmsimage.py
  shapefit.py
  shapelets.py
  sourcecounts.py
  spectralindex.py
  statusbar.py
  tc.py
  threshold.py
  wavelet_atrous.py
  DESTINATION lofar/bdsm)
