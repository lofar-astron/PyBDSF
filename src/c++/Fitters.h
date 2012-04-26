
#ifndef _FITTERS_H_INCLUDED
#define _FITTERS_H_INCLUDED

#include "MGFunction.h"

/*!
  \file Fitters.h
  
  \ingroup pybdsm
  
  \author Oleksandr Usov

  \date 30/10/2007
*/

bool lmder_fit(MGFunction &fcn, bool final, int verbose);
bool dn2g_fit(MGFunction &fcn, bool final, int verbose);
bool dnsg_fit(MGFunction &fcn, bool final, int verbose);

#endif // _FITTERS_H_INCLUDED
