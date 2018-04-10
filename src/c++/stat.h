#ifndef _CBDSM_STAT_H_INCLUDED_
#define _CBDSM_STAT_H_INCLUDED_

#include <boost/python.hpp>
#include <numpy/arrayobject.h>

#if BOOST_VERSION < 106500
  typedef typename boost::python::pyndarray pyndarray;
  //namespace np = boost::python::numeric;
#else
  #include <boost/python/numpy.hpp>
  typedef typename boost::python::numpy::ndarray pyndarray;
  namespace np = boost::python::numpy;
#endif
/*!
  \file stat.h
  
  \ingroup pybdsm
  
  \author Oleksandr Usov

  Clipped RMS and mean value calculation for numpy array.
*/

boost::python::object bstat (pyndarray arr,
			     boost::python::object mask,
			     double kappa);

#endif // _CBDSM_STAT_H_INCLUDED_
