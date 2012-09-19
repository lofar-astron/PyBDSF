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
   This fitter uses DNSG algorithm from the (free subset of) PORT3
   library (www.netlib.org/port).

   DNSG is an separable least-squares solver which calls DN2G
   algorithm to solve for non-linear parameters and then uses
   linear least-squares solver to compute the linear parameters.
**/


// this is prototype for a fortran routines
extern "C"
void dnsg_ (int &n, int &p, int &l,
	    double *alf, double *c, double *Y, void *F, void *J,
	    int *inc, int &iinc, int *iv, int &liv, int &lv, double *v,
	    int *uiparm, double *urparm, void *ufparm);

extern "C"
void divset_(const int &alg, int *iv, int &liv, int &lv, double *v);


// user functions
// FIXME: these should have been declared "extern "C" static ...", but gcc 4.2.2 rejects such declarations
static
void dnsg_f(int &n, int &p, int &l, double *alf, int &nf, double *phi,
	    void *uiparm, void *urparm, void *ufparm);

static
void dnsg_df(int &n, int &p, int &l, double *alf, int &nf, double *der,
	     void *uiparm, void *urparm, void *ufparm);


// dnsg driver
bool dnsg_fit(MGFunction &fcn, bool final, int verbose)
{
  int dsize = fcn.data_size();
  int lpar  = fcn.gaul_size();
  int nlpar = fcn.parameters_size() - lpar;

  // working variables
  int n = dsize, p = nlpar, l = lpar;
  int iinc = l+1, liv = 115 + p + l + 2*p;
  int lv = 105 + n*(l+p+3) + (l+p)*(n+l+p+1) + l*(l+3)/2 + p*(2*p+17);
  vector<double> alf(p), c(l), y(n), v(lv);
  vector<int> inc((l+1)*p), iv(liv);

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

  fcn.get_nlin_parameters(&alf[0]);
  fcn.data(&y[0]);

  // fill in incidence matrix
  int pi = 0;
  for (int li = 0; li < l; ++li) {
    for (int i = 0; i < fcn.gaussian_size(li) - 1; ++i) {
      inc[li + pi * (l+1)] = 1;
      ++ pi;
    }
  }

  dnsg_(n, p, l,
	&alf[0], &c[0], &y[0], (void *)dnsg_f, (void *)dnsg_df,
	&inc[0], iinc, &iv[0], liv, lv, &v[0],
	(int *)0, (double *)0, (void *)&fcn);

  fcn.set_nlin_parameters(&alf[0]);
  fcn.set_lin_parameters(&c[0]);

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
	 << "  DNSG"
	 << endl;
  }

  return converged;
}

// user-supplied functions
static void dnsg_f(int &n, int &p, int &l, double *alf, int &nf, double *phi,
		   void *uiparm, void *urparm, void *ufparm)
{
  (void)nf;
  (void)uiparm;
  (void)urparm;

  MGFunction *fcn = (MGFunction *)ufparm;

  assert(n == fcn->data_size());
  assert(p == fcn->parameters_size() - fcn->gaul_size());
  assert(l == fcn->gaul_size());

  fcn->set_nlin_parameters(alf);
  fcn->fcn_partial_value(phi);
}

static void dnsg_df(int &n, int &p, int &l, double *alf, int &nf, double *der,
		    void *uiparm, void *urparm, void *ufparm)
{
  (void)nf;
  (void)uiparm;
  (void)urparm;

  MGFunction *fcn = (MGFunction *)ufparm;

  assert(n == fcn->data_size());
  assert(p == fcn->parameters_size() - fcn->gaul_size());
  assert(l == fcn->gaul_size());

  fcn->set_nlin_parameters(alf);
  fcn->fcn_partial_gradient(der);
}

