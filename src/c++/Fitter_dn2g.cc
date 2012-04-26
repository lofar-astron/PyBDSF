/*!
  \file Fitter_dn2g.cc
  
  \ingroup pybdsm
  
  \author Oleksandr Usov
  
  \date 30/10/2007
*/

#include "Fitters.h"
#include <iostream>

using namespace std;

/**
   This fitter uses DN2G algorithm from the (free subset of) PORT3 
   library (www.netlib.org/port).

   DN2G is an implementation of the NL2SOL -- adaptive nonlinear least squares
   algorithm; also known as TOMS/573.
**/


// this is prototype for a fortran routines
extern "C"
void dn2g_(int &n, int &p, double *x, void *F, void *J,
	   int *iv, int &liv, int &lv, double *v,
	   int *uiparm, double *urparm, void *ufparm);

extern "C"
void divset_(const int &alg, int *iv, int &liv, int &lv, double *v);


// user functions
// FIXME: these should have been declared "extern "C" static ...", but gcc 4.2.2 rejects such declarations
static
void dn2g_f(int &n, int &p, double *x, int &nf, double *F,
	    void *uiparm, void *urparm, void *ufparm);

static
void dn2g_df(int &n, int &p, double *x, int &nf, double *J,
	     void *uiparm, void *urparm, void *ufparm);


// dn2g driver
bool dn2g_fit(MGFunction &fcn, bool final, int verbose)
{
  int dsize = fcn.data_size();
  int npar = fcn.parameters_size();

  // working variables
  int n = dsize, p = npar, liv = 82+p, lv = 105 + p*(n+2*p+17) + 2*n;
  vector<double> x(p), v(lv);
  vector<int> iv(liv);

  // set run-time parameters
  divset_(1, &iv[0], liv, lv, &v[0]);

  iv[16] = 1000;  // MXFCAL - maximal number of function calls
  iv[17] = 1000;  // MXITER - maximal number of iterations
  if (final)
    v[32] = 1e-8; // XCTOL  - x-convergence tolerance
  else
    v[32] = 1e-4; // XCTOL  - x-convergence tolerance

  verbose = (verbose < 0) ? 0 : verbose;

  switch (verbose) {
  case 0: // silence it completely
  case 1:
    iv[20] = 0; // PRUNIT
    break;

  case 2: // print out some results
    iv[13] = 0;  // COVPRT - covariance & regression diagnostic (0 ... 3)
    iv[18] = 1;  // OUTLEV - iteration summary (0, 1, ...)
    iv[19] = 1;  // PARPRT - non-default parameters (0, 1)
    iv[21] = 1;  // SOLPRT - solution (x, F, J) (0, 1)
    iv[22] = 1;  // STATPR - summary statistics (-1, 0, 1)
    iv[23] = 0;  // X0PRT  - input X (0, 1)
    break;
  }

  iv[56] = 0; // RDREQ -- do not calculate covariance & regression diagnostic arrays

  fcn.get_parameters(&x[0]);

  dn2g_(n, p, &x[0], (void *)dn2g_f, (void *)dn2g_df,
	&iv[0], liv, lv, &v[0],
	(int *)0, (double *)0, (void *)&fcn);

  fcn.set_parameters(&x[0]);

  // extract errors information
  if (final) {
    // TODO
  }

  // check convergence and (possibly) print out status line
  int info = iv[0];
  bool converged = (info > 2) && (info < 7);

  if (verbose) {
    int nfev = iv[5];  // NFCALL - number of function evaluations
    int njev = iv[29]; // NGCALL - number of gradient evaluations
    double chi2 = fcn.chi2();
    cout << "status: " << converged
	 << "  code: " << info
	 << "  Fev/Jev: " << nfev << "/" << njev
	 << "  chi2(/dp): " << chi2 << "(" << chi2/dsize << ")"
	 << "  DN2G"
	 << endl;
  }

  return converged;
}

// user-supplied functions
static void dn2g_f(int &n, int &p, double *x, int &nf, double *F,
		   void *uiparm, void *urparm, void *ufparm)
{
  MGFunction *fcn = (MGFunction *)ufparm;

  assert(n == fcn->data_size());
  assert(p == fcn->parameters_size());

  fcn->set_parameters(x);
  fcn->fcn_diff(F);
}

static void dn2g_df(int &n, int &p, double *x, int &nf, double *J,
		    void *uiparm, void *urparm, void *ufparm)
{
  MGFunction *fcn = (MGFunction *)ufparm;

  assert(n == fcn->data_size());
  assert(p == fcn->parameters_size());

  fcn->set_parameters(x);
  fcn->fcn_diff_transposed_gradient(J);
}
