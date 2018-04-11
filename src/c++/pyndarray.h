#if BOOST_VERSION < 106500
  typedef typename boost::python::numeric::array pyndarray;
  //namespace np = boost::python::numeric;
#else
  #include <boost/python/numpy.hpp>
  typedef typename boost::python::numpy::ndarray pyndarray;
  namespace np = boost::python::numpy;
#endif