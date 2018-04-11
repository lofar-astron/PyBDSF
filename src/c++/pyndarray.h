#if BOOST_VERSION < 106500
  typedef typename boost::python::numeric::array pyndarray;
#else
  #include <boost/python/numpy.hpp>
  typedef typename boost::python::numpy::ndarray pyndarray;
#endif