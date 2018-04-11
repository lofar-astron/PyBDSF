#include "Python.h"
#ifndef _MGFUNCTION_H_INCLUDED
#define _MGFUNCTION_H_INCLUDED

#include <vector>
#include <utility>
#include <boost/python.hpp>
#include <pyndarray.h>

using namespace boost::python;

/*!
  \class MGFunction
  
  \ingroup pybdsm
  
  \brief Multi-Gaussian function.

  \author Oleksandr Usov

  \date 15/10/2007
  
  This class allows you to manage multi-gaussian function and implements
  all math required to use it for fitting.

  In order to improve fitting performance a number of tricks are done.
  MGFunction maintains internal caches of unmasked data points (mm_data)
  and partially evaluated gaussians (mm_fcn). 

  End-user interface consists of functions to add (py_add_gaussian),
  remove (py_remove_gaussian), retrieve (py_get_gaussian) single gaussians
  and a functions to access all parameters at once (py_get/set_parameters).

  Few more support routines are present too (py_find_peak, py_reset, etc).


  Internal interface for fitter routines provides a number of low-level
  functions to evaluate gaussians (and their derivatives) in a number of
  ways (fcn_value, fcn_diff, fcn_gradient, etc).

  An important note -- current implementation isn't thread-safe, as caches
  are shared between all instances. One of the possibilities change it is
  to define caches thread-local, but this will require special care for the
  cleanup to prevent memory leaks.
*/

class MGFunction
{
 public:
  MGFunction(pyndarray data, pyndarray mask, double weight);
  ~MGFunction();
  
  ////////////////////////////////
  // High-level Python interface
  ////////////////////////////////
  enum Gtype {
    G_Amplitude_Only = 1,
    G_Reduced_Gaussian = 3,
    G_Gaussian = 6,
  };
  
  void py_reset();
  void py_add_gaussian(Gtype type, object parameters);
  void py_remove_gaussian(int idx);
  tuple py_get_gaussian(int idx);
  void py_set_gaussian(int idx, object parameters);
  list py_get_parameters();
  void py_set_parameters(object parameters);
  list py_get_errors();
  tuple py_find_peak();
  static void register_class();
  
  ////////////////////////////////
  // Low-level interface for fitting routines
  ////////////////////////////////
  int data_size() const { return m_ndata; }
  int gaul_size() const { return m_gaul.size(); }
  int parameters_size() const { return m_npar; }
  int gaussian_size(unsigned idx) const { assert(idx < m_gaul.size()); return m_gaul[idx]; }
  
  void get_parameters(double *buf) const;
  void set_parameters(const double *buf);
  void get_nlin_parameters(double *buf) const;
  void set_nlin_parameters(const double *buf);
  void set_lin_parameters(const double *buf);
  
  /*! data (unmasked pixels only) */
  void data(double *buf) const;
  /*! value of the function (unmasked pixels only) */
  void fcn_value(double *buf) const;
  /*! data - value (unmasked pixels only) */
  void fcn_diff(double *buf) const;
  /*! values of single (unscaled) gaussians */
  void fcn_partial_value(double *buf) const;
  /*! fcn_value gradient */
  void fcn_gradient(double *buf) const;
  /*! fcn_diff gradient */
  void fcn_diff_gradient(double *buf) const;
  /*! fcn_value gradient (transposed) */
  void fcn_transposed_gradient(double *buf) const;
  /*! fcn_diff gradient (transposed) */
  void fcn_diff_transposed_gradient(double *buf) const;
  /*! fcn_partial_value gradient */
  void fcn_partial_gradient(double *buf) const;
  /*! calculate chi^2 of the residual image */
  double chi2() const;

protected:
  
  /*! gaussians types (lengths) */
  std::vector<int> m_gaul;
  /*! parameters of gaussians */
  std::vector<std::vector<double> > m_parameters;
  /*! error bars */
  std::vector<std::vector<double> > m_errors;
  
  /*! weight for chi^2 calculation */
  double m_weight;
  /*! number of fitted parameters */
  unsigned m_npar;
  /*! number of fitted (unmasked) datapoints */
  unsigned m_ndata;
  /*! Data array */
  pyndarray m_data;
  /*! Mask array */
  pyndarray m_mask;
  
 private:
  /*! prevent copying of the MGFunction objects */
  MGFunction(MGFunction const &);
  
  /// these are used to cache intermediate calculations
  template<class T>
    void __update_dcache() const;
  void _update_dcache() const; /// update cached data
  void _update_fcache() const; /// update cached function
  unsigned long _cksum() const; /// cksum of m_parameters


  typedef struct {
    int x1, x2;
    double d;
  } dcache_t;

  typedef struct {
    double sn, cs, f1, f2, val;
  } fcache_t;

  typedef std::vector<dcache_t>::iterator dcache_it;
  typedef std::vector<fcache_t>::iterator fcache_it;

  /*! Data cache */
  static std::vector<dcache_t> mm_data;
  /*! cache for function/gradient evaluations */
  static std::vector<fcache_t> mm_fcn;

  // these are used to verify whether cached values are up-to-date
  static void *mm_obj;
  static unsigned long mm_cksum;
};

#endif // _MGFUNCTION_H_INCLUDED
