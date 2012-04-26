/*!
  \file cbdsm_main.cc
  
  \ingroup pybdsm
  
  \author Oleksandr Usov
*/

#define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle

#include "stat.h"
#include "MGFunction.h"
#include "Fitters.h"
#include <num_util/num_util.h>

using namespace boost::python;

BOOST_PYTHON_MODULE(_cbdsm)
{
  import_array();
  numeric::array::set_module_and_type("numpy", "ndarray");

  scope().attr("__doc__") = 
    "A collection of optimized C & Fortran routines for pybdsm";

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
