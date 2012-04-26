/*!
  \file Fitter_lmder.cc
  
  \ingroup pybdsm
  
  \author Oleksandr Usov
  
  \date 30/10/2007
*/

#include "Fitters.h"
#include <iostream>

using namespace std;

/**
   This fitter uses slightly modified version of the Levenberg-Marquardt
   algorithm lmder from netlib's MINPACK-1 (www.netlib.org/minpack).

   The modification of the original lmder.f was needed to pass extra
   arguments into user function.
**/


// this is prototype for a fortran minimization routine
extern "C" 
void lmder_(void *fcn, int &m, int &n, double *x, double *F, double *J, int &ldfjac,
	    double &ftol, double &xtol, double &gtol, int &maxfev,
	    double *diag, int &mode, double &factor,
	    int &nprint, int &info, int &nfev, int &njev,
	    int *ipvt, double *qtf, double *wa1, double *wa2, double *wa3, double *wa4,
	    void *userpar);

// user function
// FIXME: these should have been declared "extern "C" static ...", but gcc 4.2.2 rejects such declarations
static
void lmder_fcn(int &m, int &n, double *x, double *F, double *J, int &ldfjac, 
	       int &iflag, void *userpar);


// lmder driver
bool lmder_fit(MGFunction &fcn, bool final, int verbose)
{
  int dsize = fcn.data_size();
  int npar = fcn.parameters_size();

  // working variables
  int m = dsize, n = npar, ldfjac = m, maxfev = 200,
    mode = 1, nprint = 0, info, nfev, njev;
  double ftol, xtol, gtol, factor = 10;
  vector<double> x(n), F(m), J(ldfjac * n), diag(n),
    qtf(n), wa1(n), wa2(n), wa3(n), wa4(m);
  vector<int> ipvt(n);

  // set run-time parameters
  gtol = 0;
  if (final)
    ftol = xtol = 1e-8;
  else
    ftol = xtol = 1e-4;

  fcn.get_parameters(&x[0]);

  lmder_((void *)lmder_fcn, m, n, &x[0], &F[0], &J[0], ldfjac,
	 ftol, xtol, gtol, maxfev,
	 &diag[0], mode, factor,
	 nprint, info, nfev, njev,
	 &ipvt[0], &qtf[0], &wa1[0], &wa2[0], &wa3[0], &wa4[0],
	 (void *)&fcn);

  fcn.set_parameters(&x[0]);

  // extract errors information
  if (final) {
    // TODO
  }

  // check convergence and (possibly) print out status line
  bool converged = (info > 0) && (info < 4);

  if (verbose) {
    double chi2 = fcn.chi2();
    cout << "status: " << converged
	 << "  code: " << info
	 << "  Fev/Jev: " << nfev << "/" << njev
	 << "  chi2(/dp): " << chi2 << "(" << chi2/dsize << ")"
	 << "  LMDER"
	 << endl;
  }

  return converged;
}

// user-supplied function
static void lmder_fcn(int &m, int &n, double *x, double *F, double *J, int &ldfjac, 
		       int &iflag, void *userpar)
{
  MGFunction *fcn = (MGFunction *)userpar;

  assert(m == fcn->data_size());
  assert(n == fcn->parameters_size());

  fcn->set_parameters(x);

  switch (iflag) {
  case 1:
    return fcn->fcn_diff(F);

  case 2:
    return fcn->fcn_diff_transposed_gradient(J);

  default:
    cerr << 
      "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
      " LMDER C-wrapper\n"
      " unexpected value of iflag\n"
      "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    abort();
  }
}
