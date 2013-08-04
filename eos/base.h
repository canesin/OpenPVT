//  OpenPROP - Open source properties library
//	Fabio Cesar Canesin <fabio.canesin@gmail.com>
//	MIT 2013 license
//	base.h - Basic thermodynamic relations for OpenPROP

#ifndef BASE_H
#define BASE_H

#include <math.h>
#include <limits>
#include "../utils/utils.h"

namespace openprop
{

// W. R. Salzman, Critical Constants of the van der Waals Gas,
// Van der Waals and ideal gas EoS to provide guess for PNRB
// @return pressure from temperature, and density
REAL P_guess(const fluid& F, const REAL& T, const REAL& rho){
  //Compute reduced state
  REAL tau = T/F.Tc
  REAL delta = rho/F.rhoc
  return  delta*(-3*(delta*delta) + 9*delta - 8*tau)/(delta - 3)
}
// @return temperature from pressure, and density
REAL T_guess(const fluid& F, const REAL& P, const REAL& rho){
  //Compute reduced state
  REAL p = P/F.pc
  REAL delta = rho/F.rhoc
  return  -(delta -3)*(p + 3*(delta*delta))/(8*delta + std::numeric_limits<REAL>::min())
}
// @return density from temperature, and pressure
REAL Rho_guess(const fluid& F, const REAL& T, const REAL& P){
  //Compute reduced state
  return  R * T/(P + std::numeric_limits<REAL>::min())
}

} /* openprop */
#endif /* end of include guard: BASE_H */