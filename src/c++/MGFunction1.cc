/*!
  \file MGFunction.cc
  
  \ingroup pybdsm
  
  \author Oleksandr Usov
  
  \author 15/10/2007


Python interface.
Most routines here are pretty straighforward, as they just wrap/unwrap
parameters coming to/from python side.
*/

#define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle
#define NO_IMPORT_ARRAY

#include "boost_python.h"
#include "MGFunction.h"

#include <num_util/num_util.h>
#include <cfloat>

using namespace std;
namespace n = num_util;


//
// Constructor -- check data types/shapes and store them
//
MGFunction::MGFunction(numeric::array data, numeric::array mask, double weight)
  :  m_weight(weight), m_npar(0), m_data(data), m_mask(mask)
{
  py_assert(n::rank(data) == 2 && n::rank(mask) == 2,
            PyExc_ValueError, "Data and mask should be rank-2 arrays");
  py_assert(n::shape(data) == n::shape(mask),
	    PyExc_ValueError, "Shape of mask doesn't matches that of data");
  py_assert(n::type(mask) == NPY_BOOL,
            PyExc_TypeError, "Incorrect mask datatype");

  // to calculate m_ndata we subtract masked pixels
  PyObject *sum = PyArray_Sum((PyArrayObject *)mask.ptr(), NPY_MAXDIMS, NPY_INT, NULL);
  m_ndata = n::size(data) - PyInt_AsLong(sum);
  Py_DECREF(sum);
}

//
// Destructor -- essentially do-nothing thing
//
MGFunction::~MGFunction()
{
  // enforce data cache reset
  // this is needed if new MGFunction object is allocated 
  // at the same spot in memory as previous one
  if (mm_obj == this)
    mm_obj = 0;
}

//
// Clear MGFunction, dropping all gaussians
//
void MGFunction::py_reset()
{
  m_npar = 0;
  m_gaul.clear();
  m_parameters.clear();
  m_errors.clear();
}

//
// Add a gaussian of a specific kind
//
void MGFunction::py_add_gaussian(Gtype type, object parameters)
{
  py_assert(len(parameters) == 6,
	    PyExc_ValueError, "Wrong number of parameters for gaussian");

  vector<double> t(6);
  for (int i = 0; i < 6; ++i)
    t[i] = extract<double>(parameters[i]);

  m_npar += int(type);
  m_gaul.push_back(int(type));
  m_parameters.push_back(t);
  m_errors.push_back(vector<double>(6));
}

//
// Remove gaussian by index
//
void MGFunction::py_remove_gaussian(int idx)
{
  if (idx < 0)
    idx += m_gaul.size();

  py_assert(idx >= 0 && idx < (int)m_gaul.size(),
	    PyExc_IndexError, "Incorrect index");

  m_npar -= m_gaul[idx];
  m_gaul.erase(m_gaul.begin() + idx);
  m_parameters.erase(m_parameters.begin() + idx);
  m_errors.erase(m_errors.begin() + idx);
}

//
// Get gaussian parameters by index
//
tuple MGFunction::py_get_gaussian(int idx)
{
  if (idx < 0)
    idx += m_gaul.size();

  py_assert(idx >= 0 && idx < (int)m_gaul.size(),
	    PyExc_IndexError, "Incorrect index");

  vector<double> &p = m_parameters[idx];

  return make_tuple(p[0], p[1], p[2], p[3], p[4], p[5]);
}

//
// Set gaussian parameters by index
//
void MGFunction::py_set_gaussian(int idx, object parameters)
{
  if (idx < 0)
    idx += m_gaul.size();

  py_assert(idx >= 0 && idx < (int)m_gaul.size(),
	    PyExc_IndexError, "Incorrect index");
  py_assert(len(parameters) == 6,
	    PyExc_ValueError, "Wrong number of parameters for gaussian");

  for (int i = 0; i < 6; ++i)
    m_parameters[idx][i] = extract<double>(parameters[i]);
}

//
// Get all gaussian parameters as a list of tuples
//
list MGFunction::py_get_parameters()
{
  list res;

  for (unsigned i = 0; i < m_gaul.size(); ++i)
    res.append(py_get_gaussian(i));


  return res;
}

//
// Set all gaussian parameters from a list of tuples
//
void MGFunction::py_set_parameters(object parameters)
{
  py_assert(len(parameters) == int(m_gaul.size()),
	    PyExc_ValueError, "Wrong number of gaussians");

  for (unsigned i = 0; i < m_parameters.size(); ++i)
    py_set_gaussian(i, parameters[i]);
}

//
// Errors -- not really implemented now
//
list MGFunction::py_get_errors()
{
  list res;

  for (unsigned i = 0; i < m_gaul.size(); ++i) {
    vector<double> &e = m_errors[i];
    res.append(make_tuple(e[0], e[1], e[2], e[3], e[4], e[5]));
  }

  return res;
}

//
// Find highest peak in the data-MGFunction residual
//
tuple MGFunction::py_find_peak()
{
  vector<double> buf(data_size());
  fcn_diff(&buf.front());

  double peak = buf[0];
  unsigned pidx = 0;

  for (unsigned i = 0; i < buf.size(); ++i)
    if (buf[i] > peak) {
      peak = buf[i];
      pidx = i;
    }

  int x1 = mm_data[pidx].x1;
  int x2 = mm_data[pidx].x2;

  return make_tuple(peak, make_tuple(x1, x2));
}


//
// Register MGFunction class in python interpreter
//
void MGFunction::register_class()
{
  /*
    FIXME: I don't really understand why this 'using' statement should be here, but
    my compiler (apple gcc 4.0.1) can't find arg without it.
   */
  using boost::python::arg;

  enum_<Gtype>("Gtype")
    .value("g1", G_Amplitude_Only)
    .value("g3", G_Reduced_Gaussian)
    .value("g6", G_Gaussian)
    ;

  class_<MGFunction,
    boost::noncopyable>("MGFunction",
			"Multi-Gaussian function.\n\n"
			"This class allows you to manage multi-gaussian function\n"
			"and implements all math required to use it for fitting.\n\n"
			"NEVER EVER USE IT IN MULTITHREADED SOFTWARE WITHOUT APPROPRIATE LOCKING\n"
			"IT'S INTERNAL CACHES ARE NOT THREAD-SAFE\n\n",
			init<numeric::array, numeric::array, 
			double>((arg("data"), "mask", arg("weight") = 1.)))

    .def("__len__", &MGFunction::gaul_size, 
	 "total number of gaussians")

    .def("add_gaussian", &MGFunction::py_add_gaussian, (arg("type"), "parameters"),
	 "add gaussian to the set")

    .def("remove_gaussian", &MGFunction::py_remove_gaussian, arg("idx"),
	 "remove given gaussian")

    .def("__delitem__", &MGFunction::py_remove_gaussian, arg("idx"),
	 "remove given gaussian")

    .def("__getitem__", &MGFunction::py_get_gaussian, arg("idx"),
	 "get parameters of the given gaussian")

    .def("__setitem__", &MGFunction::py_set_gaussian, (arg("idx"), "parameters"),
	 "set parameters of the given gaussian")

    .def("find_peak", &MGFunction::py_find_peak,
	 "find highest peak in the residual image\n"
	 "returns (value, (x1, x2))")

    .def("fitted_parameters", &MGFunction::parameters_size,
	 "total number of fitted parameters")

    .def("reset", &MGFunction::py_reset,
	 "reset MGFunction by forgetting all gaussians")

    ADD_PROPERTY2("parameters", &MGFunction::py_get_parameters, &MGFunction::py_set_parameters,
		  "get/set parameters of all gaussians together")

    ADD_PROPERTY1("errors", &MGFunction::py_get_errors,
		  "get error bounds")
    ;
}
