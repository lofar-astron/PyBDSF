## ----------------------------------------------------------------------
## $Id$
## ----------------------------------------------------------------------

enable_language(Fortran)

set (minpack_sources
     lmder.f lmpar.f qrfac.f qrsolv.f enorm.f dpmpar.f)

add_definitions (-fPIC)

add_library (minpack STATIC ${minpack_sources})

