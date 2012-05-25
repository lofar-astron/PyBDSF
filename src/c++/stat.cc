/*!
  \file stat.cc
  
  \ingroup pybdsm
  
  \author Oleksandr Usov
*/

#define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle
#define NO_IMPORT_ARRAY

#include "boost_python.h"
#include "stat.h"

#include <num_util/num_util.h>
#include <cfloat>

using namespace boost::python;
using namespace std;
namespace n = num_util;

// check if obj is numpy scalar object and it's value is val
static bool
npybool_check(PyObject *obj, bool val)
{
  npy_bool tmp;

  if (PyArray_IsScalar(obj, Bool)) {
    PyArray_ScalarAsCtype(obj, &tmp);
    return tmp == val;
  }

  return false;
}

static inline
bool _nonzero(const unsigned &n, const int *v)
{
  bool res = true;
  for (unsigned i = 0; i < n; ++i)
    res = res && v[i];

  return res;
}

// calculate clipped mean/rms for array
template<class T>
static pair<double, double>
_stat_nd(numeric::array arr, double _mean, double _threshold)
{
  // ensure contiguous memory access by appropriately sorting indices
  vector<int> shape = n::shape(arr);
  vector<int> strides = n::strides(arr);
  const unsigned Nd = shape.size();

  for (unsigned i = 0; i < Nd; ++i)
    strides[i] /= sizeof(T);

  for (unsigned i = 0; i < Nd; ++i)
    for (unsigned j = i; j < Nd; ++j)
      if (strides[i] > strides[j]) {
	swap(shape[i], shape[j]);
	swap(strides[i], strides[j]);
      }

  // accumulators
  double sum = 0;
  double sumsq = 0;
  int N = n::size(arr);

  // data access machinery: we need adv to simplify jumping from
  // one line to the next one (offset from the end of previous line)
  T *data = (T *)n::data(arr);
  int idx[Nd];
  int adv[Nd];

  for (unsigned i = 0; i < Nd; ++i) {
    idx[i] = shape[i];
    adv[i] = strides[i];
    if (i)
      adv[i] -= shape[i-1]*strides[i-1];
  }

  while (_nonzero(Nd, idx)) {
    // inner loop over fastest index
    int _cnt = idx[0];
    int _adv = adv[0];
    while(_cnt) {
      --_cnt;
      double v = *data;
      if (fabs(v - _mean) > _threshold)
	--N;
      else {
	sum += v;
	sumsq += v*v;
      }
      data += _adv;
    }
    idx[0] = _cnt;
    
    // now a bit of magic to handle all other indices
    for (unsigned i = 1; i < Nd; ++i) {
      idx[i-1] = shape[i-1];
      data += adv[i];
      --idx[i];
      if (idx[i] != 0)
	break;
    }
  }

  double mean = sum/N;
  double dev  = sqrt((sumsq + N*mean*mean - 2*mean*sum)/(N-1));

  return make_pair(mean, dev);
}


// calculate clipped mean/rms for masked array
template<class T>
static pair<double, double>
_stat_nd_m(numeric::array arr, numeric::array mask, double _mean, double _threshold)
{
  // ensure contiguous memory access by appropriately sorting indices
  vector<int> shape = n::shape(arr);
  vector<int> strides = n::strides(arr);
  vector<int> mstrides = n::strides(mask);
  const unsigned Nd = shape.size();

  for (unsigned i = 0; i < Nd; ++i) {
    strides[i] /= sizeof(T);
    mstrides[i] /= sizeof(npy_bool);
  }

  for (unsigned i = 0; i < Nd; ++i)
    for (unsigned j = i; j < Nd; ++j)
      if (strides[i] > strides[j]) {
	swap(shape[i], shape[j]);
	swap(strides[i], strides[j]);
	swap(mstrides[i], mstrides[j]);
      }

  // accumulators
  double sum = 0;
  double sumsq = 0;
  int N = n::size(arr);

  // data access machinery: we need adv to simplify jumping from
  // one line to the next one (offset from the end of previous line)
  T *data = (T *)n::data(arr);
  npy_bool *mdata = (npy_bool *)n::data(mask);
  int idx[Nd];
  int adv[Nd];
  int madv[Nd];

  for (unsigned i = 0; i < Nd; ++i) {
    idx[i] = shape[i];
    adv[i] = strides[i];
    madv[i] = mstrides[i];
    if (i) {
      adv[i] -= shape[i-1]*strides[i-1];
      madv[i] -= shape[i-1]*mstrides[i-1];
    }
  }

  while (_nonzero(Nd, idx)) {
    // inner loop over fastest index
    int _cnt = idx[0];
    int _adv = adv[0];
    int _madv = madv[0];
    while(_cnt) {
      --_cnt;
      double v = *data;
      if (*mdata || fabs(v - _mean) > _threshold)
	--N;
      else {
	sum += v;
	sumsq += v*v;
      }
      data += _adv;
      mdata += _madv;
    }
    idx[0] = _cnt;
    
    // now a bit of magic to handle all other indices
    for (unsigned i = 1; i < Nd; ++i) {
      idx[i-1] = shape[i-1];
      data += adv[i];
      mdata += madv[i];
      --idx[i];
      if (idx[i] != 0)
	break;
    }
  }

  double mean = sum/N;
  double dev  = sqrt((sumsq + N*mean*mean - 2*mean*sum)/(N-1));

  return make_pair(mean, dev);
}

// dispatch calculation to the correct _stat_FOO function
template<class T>
static pair<double, double>
_stat(numeric::array arr, object mask, double _mean, double _threshold)
{
  if (mask.ptr() == Py_None || mask.ptr() == Py_False 
      || npybool_check(mask.ptr(), false))
    return _stat_nd<T>(arr, _mean, _threshold);
  else
    return _stat_nd_m<T>(arr, extract<numeric::array>(mask), _mean, _threshold);
}

template<class T>
static object _bstat(numeric::array arr, object mask, double kappa)
{
  #include <iostream>
  const int max_iter = 200;
  vector<double> mean, dev;
  pair<double, double> res(0, DBL_MAX);
  int cnt = 0;

  do {
    ++cnt;
    mean.push_back(res.first);
    dev.push_back(res.second);
    res = _stat<T>(arr, mask, mean.back(), kappa*dev.back());
  } while (res.second != dev.back() && cnt < max_iter);

 /* py_assert(res.second == dev.back(),
	    PyExc_RuntimeError, "clipped rRMS calculation does not converge"); */

  return make_tuple(mean[1], dev[1], mean.back(), dev.back(), cnt);
}

object bstat(numeric::array arr, object mask, double kappa)
{
  NPY_TYPES type = n::type(arr);

  if (PyArray_ISBYTESWAPPED(arr.ptr()))
    goto fail;

  if (mask.ptr() != Py_None && mask.ptr() != Py_False 
      && !npybool_check(mask.ptr(), false)) {
    numeric::array amask = extract<numeric::array>(mask);

    n::check_type(amask, NPY_BOOL);
    int rank = n::rank(arr);
    n::check_rank(amask, rank);
    n::check_size(amask, n::size(arr));
    // this is pointless on pc, but might matter somewhere else
    if (PyArray_ISBYTESWAPPED(amask.ptr()))
      goto fail;
  }

  switch (type) {
  case NPY_DOUBLE:
    return _bstat<npy_double>(arr, mask, kappa);
  case NPY_FLOAT:
    return _bstat<npy_float>(arr, mask, kappa);
  default:
    goto fail;
  }

 fail:
  py_assert(false, 
	    PyExc_RuntimeError, "bstat dispatch failed: not implemented for this datatype/layout");
  return tuple(); // this is fake return-statement to silence the compiler
}
