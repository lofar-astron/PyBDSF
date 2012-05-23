#ifndef _AUX_H_INCLUDED
#define _AUX_H_INCLUDED

/*!
	\file boost_python.h
	
	\ingroup pybdsm

	\brief Miscellaneous usefull routines
*/

#include <boost/version.hpp>
#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>

#if BOOST_VERSION > 103200
#define ADD_PROPERTY1(name, get, doc) .add_property(name, get, doc)
#define ADD_PROPERTY2(name, get, set, doc) .add_property(name, get, set, doc)
#else
#define ADD_PROPERTY1(name, get, doc) .add_property(name, get)
#define ADD_PROPERTY2(name, get, set, doc) .add_property(name, get, set)
#endif

inline void py_assert(bool cond, PyObject *exc, const char *msg)
{
  if(!cond) {
    PyErr_SetString(exc, msg);
    throw boost::python::error_already_set();
  }
}

#endif // _AUX_H_INCLUDED
