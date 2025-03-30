// Copyright 2006  Phil Austin (http://www.eos.ubc.ca/personal/paustin)
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle
#define NO_IMPORT_ARRAY
#include "num_util.h"

// namespace { const char* rcsid = "$Id$"; }

using namespace boost::python;

namespace num_util{

  //specializations for use by makeNum


  template <>
  NPY_TYPES getEnum<unsigned char>(void)
  {
    return NPY_UBYTE;
  }


  template <>
  NPY_TYPES getEnum<signed char>(void)
  {
    return NPY_BYTE;
  }

  template <>
  NPY_TYPES getEnum<short>(void)
  {
    return NPY_SHORT;
  }

  template <>
  NPY_TYPES getEnum<unsigned short>(void)
  {
    return NPY_USHORT;
  }


  template <>
  NPY_TYPES getEnum<unsigned int>(void)
  {
    return NPY_UINT;
  }

  template <>
  NPY_TYPES getEnum<int>(void)
  {
    return NPY_INT;
  }

  template <>
  NPY_TYPES getEnum<long>(void)
  {
    return NPY_LONG;
  }

  template <>
  NPY_TYPES getEnum<unsigned long>(void)
  {
    return NPY_ULONG;
  }


  template <>
  NPY_TYPES getEnum<long long>(void)
  {
    return NPY_LONGLONG;
  }

  template <>
  NPY_TYPES getEnum<unsigned long long>(void)
  {
    return NPY_ULONGLONG;
  }

  template <>
  NPY_TYPES getEnum<float>(void)
  {
    return NPY_FLOAT;
  }

  template <>
  NPY_TYPES getEnum<double>(void)
  {
    return NPY_DOUBLE;
  }

  template <>
  NPY_TYPES getEnum<long double>(void)
  {
    return NPY_LONGDOUBLE;
  }

  template <>
  NPY_TYPES getEnum<std::complex<float> >(void)
  {
    return NPY_CFLOAT;
  }


  template <>
  NPY_TYPES getEnum<std::complex<double> >(void)
  {
    return NPY_CDOUBLE;
  }

  template <>
  NPY_TYPES getEnum<std::complex<long double> >(void)
  {
    return NPY_CLONGDOUBLE;
  }


typedef KindStringMap::value_type  KindStringMapEntry;
KindStringMapEntry kindStringMapEntries[] =
  {
    KindStringMapEntry(NPY_UBYTE,  "NPY_UBYTE"),
    KindStringMapEntry(NPY_BYTE,   "NPY_BYTE"),
    KindStringMapEntry(NPY_SHORT,  "NPY_SHORT"),
    KindStringMapEntry(NPY_INT,    "NPY_INT"),
    KindStringMapEntry(NPY_LONG,   "NPY_LONG"),
    KindStringMapEntry(NPY_FLOAT,  "NPY_FLOAT"),
    KindStringMapEntry(NPY_DOUBLE, "NPY_DOUBLE"),
    KindStringMapEntry(NPY_CFLOAT, "NPY_CFLOAT"),
    KindStringMapEntry(NPY_CDOUBLE,"NPY_CDOUBLE"),
    KindStringMapEntry(NPY_OBJECT, "NPY_OBJECT"),
    KindStringMapEntry(NPY_NOTYPE ,"NPY_NOTYPE")
  };

typedef KindCharMap::value_type  KindCharMapEntry;
KindCharMapEntry kindCharMapEntries[] =
  {
    KindCharMapEntry(NPY_UBYTE,  'B'),
    KindCharMapEntry(NPY_BYTE,   'b'),
    KindCharMapEntry(NPY_SHORT,  'h'),
    KindCharMapEntry(NPY_INT,    'i'),
    KindCharMapEntry(NPY_LONG,   'l'),
    KindCharMapEntry(NPY_FLOAT,  'f'),
    KindCharMapEntry(NPY_DOUBLE, 'd'),
    KindCharMapEntry(NPY_CFLOAT, 'F'),
    KindCharMapEntry(NPY_CDOUBLE,'D'),
    KindCharMapEntry(NPY_OBJECT, 'O')
  };

typedef KindTypeMap::value_type  KindTypeMapEntry;
KindTypeMapEntry kindTypeMapEntries[] =
  {
    KindTypeMapEntry('B',NPY_UBYTE),
    KindTypeMapEntry('b',NPY_BYTE),
    KindTypeMapEntry('h',NPY_SHORT),
    KindTypeMapEntry('i',NPY_INT),
    KindTypeMapEntry('l',NPY_LONG),
    KindTypeMapEntry('f',NPY_FLOAT),
    KindTypeMapEntry('d',NPY_DOUBLE),
    KindTypeMapEntry('F',NPY_CFLOAT),
    KindTypeMapEntry('D',NPY_CDOUBLE),
    KindTypeMapEntry('O',NPY_OBJECT)
  };

int numStringEntries = sizeof(kindStringMapEntries)/sizeof(KindStringMapEntry);
int numCharEntries = sizeof(kindCharMapEntries)/sizeof(KindCharMapEntry);
int numTypeEntries = sizeof(kindTypeMapEntries)/sizeof(KindTypeMapEntry);

static KindStringMap kindstrings(kindStringMapEntries,
                                   kindStringMapEntries + numStringEntries);

static KindCharMap kindchars(kindCharMapEntries,
                                   kindCharMapEntries + numCharEntries);

static KindTypeMap kindtypes(kindTypeMapEntries,
                                   kindTypeMapEntries + numTypeEntries);

//Create a numarray referencing Python sequence object
pyndarray makeNum(object x){
  if (!PySequence_Check(x.ptr())){
    PyErr_SetString(PyExc_ValueError, "expected a sequence");
    throw_error_already_set();
  }
  object obj(handle<>
	     (PyArray_ContiguousFromObject(x.ptr(),NPY_NOTYPE,0,0)));
  check_PyArrayElementType(obj);
  return extract<pyndarray>(obj);
}

//Create a one-dimensional Numeric array of length n and Numeric type t
pyndarray makeNum(ssize_t n, NPY_TYPES t=NPY_DOUBLE){
  object obj(handle<>(PyArray_SimpleNew(1, &n, t)));
  return extract<pyndarray>(obj);
}

//Create a Numeric array with dimensions dimens and Numeric type t
pyndarray makeNum(std::vector<ssize_t> dimens,
		       NPY_TYPES t=NPY_DOUBLE){
  object obj(handle<>(PyArray_SimpleNew(dimens.size(), &dimens[0], t)));
  return extract<pyndarray>(obj);
}

pyndarray makeNum(const pyndarray& arr){
  //Returns a reference of arr by calling pyndarray copy constructor.
  //The copy constructor increases arr's reference count.
  return pyndarray(arr);
}

NPY_TYPES type(pyndarray arr){
  return NPY_TYPES(PyArray_TYPE((PyArrayObject*)arr.ptr()));
}

void check_type(pyndarray arr,
		NPY_TYPES expected_type){
  NPY_TYPES actual_type = type(arr);
  if (actual_type != expected_type) {
    std::ostringstream stream;
    stream << "expected Numeric type " << kindstrings[expected_type]
	   << ", found Numeric type " << kindstrings[actual_type] << std::ends;
    PyErr_SetString(PyExc_TypeError, stream.str().c_str());
    throw_error_already_set();
  }
  return;
}

//Return the number of dimensions
int rank(pyndarray arr){
  //std::cout << "inside rank" << std::endl;
  if(!PyArray_Check(arr.ptr())){
    PyErr_SetString(PyExc_ValueError, "expected a PyArrayObject");
    throw_error_already_set();
  }
  return PyArray_NDIM((PyArrayObject*)arr.ptr());
}

void check_rank(pyndarray arr, int expected_rank){
  int actual_rank = rank(arr);
  if (actual_rank != expected_rank) {
    std::ostringstream stream;
    stream << "expected rank " << expected_rank
	   << ", found rank " << actual_rank << std::ends;
    PyErr_SetString(PyExc_RuntimeError, stream.str().c_str());
    throw_error_already_set();
  }
  return;
}

int size(pyndarray arr)
{
  if(!PyArray_Check(arr.ptr())){
    PyErr_SetString(PyExc_ValueError, "expected a PyArrayObject");
    throw_error_already_set();
  }
  return PyArray_Size(arr.ptr());
}

void check_size(pyndarray arr, int expected_size){
  int actual_size = size(arr);
  if (actual_size != expected_size) {
    std::ostringstream stream;
    stream << "expected size " << expected_size
	   << ", found size " << actual_size << std::ends;
    PyErr_SetString(PyExc_RuntimeError, stream.str().c_str());
    throw_error_already_set();
  }
  return;
}

std::vector<int> shape(pyndarray arr){
  std::vector<int> out_dims;
  if(!PyArray_Check(arr.ptr())){
    PyErr_SetString(PyExc_ValueError, "expected a PyArrayObject");
    throw_error_already_set();
  }
  npy_intp* dims_ptr = PyArray_DIMS((PyArrayObject*)arr.ptr());
  int the_rank = rank(arr);
  for (int i = 0; i < the_rank; i++){
    out_dims.push_back(*(dims_ptr + i));
  }
  return out_dims;
}

int get_dim(pyndarray arr, int dimnum){
  assert(dimnum >= 0);
  int the_rank=rank(arr);
  if(the_rank < dimnum){
    std::ostringstream stream;
    stream << "Error: asked for length of dimension ";
    stream << dimnum << " but rank of array is " << the_rank << std::ends;
    PyErr_SetString(PyExc_RuntimeError, stream.str().c_str());
    throw_error_already_set();
  }
  std::vector<int> actual_dims = shape(arr);
  return actual_dims[dimnum];
}

void check_shape(pyndarray arr, std::vector<int> expected_dims){
  std::vector<int> actual_dims = shape(arr);
  if (actual_dims != expected_dims) {
    std::ostringstream stream;
    stream << "expected dimensions " << vector_str(expected_dims)
	   << ", found dimensions " << vector_str(actual_dims) << std::ends;
    PyErr_SetString(PyExc_RuntimeError, stream.str().c_str());
    throw_error_already_set();
  }
  return;
}

void check_dim(pyndarray arr, int dimnum, int dimsize){
  std::vector<int> actual_dims = shape(arr);
  if(actual_dims[dimnum] != dimsize){
    std::ostringstream stream;
    stream << "Error: expected dimension number ";
    stream << dimnum << " to be length " << dimsize;
    stream << ", but found length " << actual_dims[dimnum]  << std::ends;
    PyErr_SetString(PyExc_RuntimeError, stream.str().c_str());
    throw_error_already_set();
  }
  return;
}

bool iscontiguous(pyndarray arr)
{
  //  return arr.iscontiguous();
  return PyArray_ISCONTIGUOUS((PyArrayObject*)arr.ptr());
}

void check_contiguous(pyndarray arr)
{
  if (!iscontiguous(arr)) {
    PyErr_SetString(PyExc_RuntimeError, "expected a contiguous array");
    throw_error_already_set();
  }
  return;
}

void* data(pyndarray arr){
  if(!PyArray_Check(arr.ptr())){
    PyErr_SetString(PyExc_ValueError, "expected a PyArrayObject");
    throw_error_already_set();
  }
  return PyArray_DATA((PyArrayObject*)arr.ptr());
}

//Copy data into the array
void copy_data(pyndarray arr, char* new_data){
  char* arr_data = (char*) data(arr);
  int nbytes = PyArray_NBYTES((PyArrayObject*)arr.ptr());
  for (int i = 0; i < nbytes; i++) {
    arr_data[i] = new_data[i];
  }
  return;
}

//Return a clone of this array
pyndarray clone(pyndarray arr){
  object obj(handle<>(PyArray_NewCopy((PyArrayObject*)arr.ptr(),NPY_CORDER)));
  return makeNum(obj);
}


//Return a clone of this array with a new type
pyndarray astype(pyndarray arr, NPY_TYPES t){
  #if BOOST_VERSION < 106500
    return (pyndarray) arr.astype(type2char(t));
  #else
    return (pyndarray) arr.astype(type2dtype(type2char(t)));
  #endif
}

std::vector<int> strides(pyndarray arr){
  std::vector<int> out_strides;
  if(!PyArray_Check(arr.ptr())){
    PyErr_SetString(PyExc_ValueError, "expected a PyArrayObject");
    throw_error_already_set();
  }
  npy_intp* strides_ptr = PyArray_STRIDES((PyArrayObject*)arr.ptr());
  int the_rank = rank(arr);
  for (int i = 0; i < the_rank; i++){
    out_strides.push_back(*(strides_ptr + i));
  }
  return out_strides;
}

int refcount(pyndarray arr){
  return Py_REFCNT(arr.ptr());
}

void check_PyArrayElementType(object newo){
  NPY_TYPES theType=NPY_TYPES(PyArray_TYPE((PyArrayObject*)newo.ptr()));
  if(theType == NPY_OBJECT){
      std::ostringstream stream;
      stream << "array elments have been cast to NPY_OBJECT, "
             << "numhandle can only accept arrays with numerical elements"
	     << std::ends;
      PyErr_SetString(PyExc_TypeError, stream.str().c_str());
      throw_error_already_set();
  }
  return;
}

std::string type2string(NPY_TYPES t_type){
  return kindstrings[t_type];
}

char type2char(NPY_TYPES t_type){
  return kindchars[t_type];
}

#if BOOST_VERSION >= 106500
boost::python::numpy::dtype type2dtype(char t){

  switch(t) {
    case 'B':
      return boost::python::numpy::dtype::get_builtin<unsigned char>();
    case 'b':
      return boost::python::numpy::dtype::get_builtin<signed char>();
    case 'h':
      return boost::python::numpy::dtype::get_builtin<short>();
    case 'i':
      return boost::python::numpy::dtype::get_builtin<int>();
    case 'l':
      return boost::python::numpy::dtype::get_builtin<long int>();
    case 'f':
      return boost::python::numpy::dtype::get_builtin<float>();
    case 'd':
      return boost::python::numpy::dtype::get_builtin<double>();
    case 'F':
      return boost::python::numpy::dtype::get_builtin<std::complex<float> >();
    case 'D':
      return boost::python::numpy::dtype::get_builtin<std::complex<double> >();
    default:
      std::cout << "Invalid character code!" << std::endl;
      break;
  }
}
#endif

NPY_TYPES char2type(char e_type){
  return kindtypes[e_type];
}

template <class T>
inline std::string vector_str(const std::vector<T>& vec)
{
  std::ostringstream stream;
  stream << "(" << vec[0];

  for(std::size_t i = 1; i < vec.size(); i++){
    stream << ", " << vec[i];
  }
  stream << ")";
  return stream.str();
}

inline void check_size_match(std::vector<int> dims, int n)
{
  int total = std::accumulate(dims.begin(),dims.end(),1,std::multiplies<int>());
  if (total != n) {
    std::ostringstream stream;
    stream << "expected array size " << n
           << ", dimensions give array size " << total << std::ends;
    PyErr_SetString(PyExc_TypeError, stream.str().c_str());
    throw_error_already_set();
  }
  return;
}


} //namespace num_util
