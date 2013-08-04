//  OpenPROP - Open source properties library
//	Fabio Cesar Canesin <fabio.canesin@gmail.com>
//	MIT 2013 license
//	eos.h - Equations of state for OpenPROP

#ifndef EOS_H
#define EOS_H

#include <math.h>
#include "../utils/utils.h"
#include "base.h"

namespace openprop
{

//Fluid Phase Equilibria 222–223 (2004) 107–118
//Lixin Sun, James F. Ely	
//Universal equation of state for engineering application:
//algorithm and application to non-polar and polar fluids
// @return Helmholtz free energy alpha
REAL SunEly(const fluid& F, const REAL& T, const REAL& rho){
  //Compute reduced state
  REAL tau = T/F.Tc
  REAL delta = rho/F.rhoc
  return  (F.SE1)*delta*pow(tau,1.5) + \
          (F.SE2)*delta*pow(tau,0.25) + \
          (F.SE3)*delta*pow(tau,1.25) + \
          (F.SE4)*pow(delta,3.0)*pow(tau,0.25) + \
          (F.SE5)*pow(delta,7.0)*pow(tau,0.875) + \
          (F.SE6)*pow(delta,2.0)*pow(tau,1.375) + \
          (F.SE7)*delta*exp(-delta) + \
          (F.SE8)*delta*pow(tau,2.375)*exp(-delta) + \
          (F.SE9)*pow(delta,2.0)*pow(tau,2.0)*exp(-delta) + \
          (F.SE10)*pow(delta,5.0)*pow(tau,2.125)*exp(-delta) + \
          (F.SE11)*delta*pow(tau,3.5)*exp(-(pow(delta,2))) + \
          (F.SE12)*delta*pow(tau,6.5)*exp(-(pow(delta,2))) + \
          (F.SE13)*(pow(delta,4))*(pow(tau,4.75))*(exp(-(pow(delta,2)))) + \
          (F.SE14)*(pow(delta,2))*(pow(tau,12.5))*(exp(-(pow(delta,3))));
}

} /* openprop */
#endif /* end of include guard: EOS_H */