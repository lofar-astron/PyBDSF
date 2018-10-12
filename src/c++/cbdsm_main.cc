/*!
  \file cbdsm_main.cc

  \ingroup pybdsf

  \author Oleksandr Usov
*/

#define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle

#include "stat.h"
#include "MGFunction.h"
#include "Fitters.h"
#include <num_util/num_util.h>

using namespace boost::python;

#if PY_MAJOR_VERSION == 2
static void wrap_import_array()
{
  import_array();
}
#else
static void * wrap_import_array()
{
  import_array();
  return NULL;
}
#endif

BOOST_PYTHON_MODULE(_cbdsm)
{
  wrap_import_array();
  #if BOOST_VERSION < 106500
  numeric::array::set_module_and_type("numpy", "ndarray");
  #else
  boost::python::numpy::initialize();
  #endif

  scope().attr("__doc__") =
    "A collection of optimized C & Fortran routines for pybdsf";

  def("bstat", &bstat, (arg("array"), arg("mask") = false, arg("kappa") = 3),
      "calculate (clipped) mean and rms of the n-dimensional (masked) image\n"
      "returns 4-tuple (mean, dev, cmean, cdev)\n");

  MGFunction::register_class();

  def("lmder_fit", &lmder_fit, (arg("fcn"), arg("final") = false, arg("verbose") = 1),
      "Fitter using the Levenberg-Marquardt algorithm LMDER from MINPACK-1");

  def("dn2g_fit", &dn2g_fit, (arg("fcn"), arg("final") = false, arg("verbose") = 1),
      "Fitter using DN2G algorithm from PORT3 library");

  def("dnsg_fit", &dnsg_fit, (arg("fcn"), arg("final") = false, arg("verbose") = 1),
      "Fitter using DNSG algorithm from PORT3 library");
}
