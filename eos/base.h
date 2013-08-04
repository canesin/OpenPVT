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
  const REAL tau = T/F.Tc
  const REAL delta = rho/F.rhoc
  return  delta*(-3*(delta*delta) + 9*delta - 8*tau)/(delta - 3)
}
// @return temperature from pressure, and density
REAL T_guess(const fluid& F, const REAL& P, const REAL& rho){
  //Compute reduced state
  const REAL p = P/F.pc
  const REAL delta = rho/F.rhoc
  return  -(delta -3)*(p + 3*(delta*delta))/(8*delta + std::numeric_limits<REAL>::min())
}
// @return density from temperature, and pressure
REAL Rho_guess(const fluid& F, const REAL& T, const REAL& P){
  //Compute reduced state
  return  R * T/(P + std::numeric_limits<REAL>::min())
}

// International Journal of Thermophysics, Vol. 24, No. 1, January 2003
// R. Span and W. Wagner
// Equations of State for Technical Applications.
// I. Simultaneously Optimized Functional
// Forms for Nonpolar and Polar Fluids
// Maxwell relations in terms of Helmholtz free energy
REAL a(const REAL& T, const REAL& rho)
{
  REAL tau = F.Tc/T;
  REAL delta = rho/F.rhoc;
  REAL R = F.R/F.MW;
  return R*T*(F.a0(tau, delta) + F.ar(tau, delta));
}

// @return P [kPa] from T [K] and rho [g/cm3]
REAL P(const REAL& T, const REAL& rho)
{
  const REAL tau = F.Tc/T;
  const REAL delta = rho/F.rhoc;
  const REAL Rf = R/F.MW;
  return Rf*T*rho*(1+delta*F.ar_d(tau, delta));
}

REAL cv(const REAL& T, const REAL& rho)
{
  const REAL tau = F.Tc/T;
  const REAL delta = rho/F.rhoc;
  const REAL Rf = R/F.MW;
  return (-pow(tau,2))*Rf*(F.a0_tt(tau, delta)+F.ar_tt(tau, delta));
}

REAL h(const REAL& T, const REAL& rho)
{
  const REAL tau = F.Tc/T;
  const REAL delta = rho/F.rhoc;
  const REAL Rf = R/F.MW;
  return Rf*T*(tau*(F.a0_t(tau, delta)+F.ar_t(tau, delta))+1+delta*F.ar_d(tau, delta));
}

REAL s(const REAL& T, const REAL& rho)
{
  REAL tau = F.Tc/T;
  REAL delta = rho/F.rhoc;
  const REAL Rf = R/F.MW;
  return Rf*(tau*(F.a0_t(tau, delta)+F.ar_t(tau, delta))-(F.ar(tau, delta)+F.a0(tau, delta)));
}

REAL u(const REAL& T, const REAL& rho)
{
  const REAL tau = F.Tc/T;
  const REAL delta = rho/F.rhoc;
  const REAL Rf = R/F.MW;
  return Rf*T*tau*(F.a0_t(tau, delta)+F.ar_t(tau, delta));
}

REAL cp(const REAL& T, const REAL& rho)
{
  const REAL tau = F.Tc/T;
  const REAL delta = rho/F.rhoc;
  const REAL Rf = R/F.MW;
  return Rf*((-pow(tau,2))*(F.a0_tt(tau, delta)+F.ar_tt(tau, delta))+ \
         pow((1+delta*F.ar_d(tau, delta)-delta*tau*F.ar_dt(tau, delta)),2)/ \
         (1+2*delta*F.ar_d(tau, delta)+pow(delta,2)*F.ar_dd(tau, delta)));
}

REAL w(const REAL& T, const REAL& rho)
{
  const REAL tau = F.Tc/T;
  const REAL delta = rho/F.rhoc;
  const REAL Rf = R/F.MW;
  return pow((Rf*T*(1+2*delta*F.ar_d(tau, delta)+(pow(delta,2))*F.ar_dd(tau, delta) - \
              pow((1+delta*F.ar_d(tau, delta)-delta*tau*F.ar_dt(tau, delta)),2)/ \
              ((pow(tau,2))*(F.a0_tt(tau, delta)+F.ar_tt(tau, delta))))),0.5);
}

REAL g(const REAL& T, const REAL& rho)
{

  return (F.a(T, rho)+F.P(T, rho)/rho);
}

REAL rho(const REAL& T, const REAL& Pressure)
{
  auto Func = [&](const REAL& rho){ return F.P(T, rho) - Pressure;};
  auto DFunc = [&](const REAL& rho){ return (Func(rho+1e-03)-Func(rho))/1e-03;};
  return F.NRB(Func, DFunc, 1e-10, 1e+10, 1e-04);
}

REAL eta0(const REAL& T)
{
  const REAL tau = T/F.Tc;
  return (F.d_0 + F.d_1*pow(tau,0.25) + F.d_2*pow(tau,0.5) + F.d_3*pow(tau,0.75))*1e-3; //Dilute gas viscosity
}

REAL eta(const REAL& T, const REAL& rho)
{
  const REAL Rf = R/F.MW;
  const REAL tau = F.Tc/T;
  const REAL delta = rho/F.rhoc;
  const REAL psi1 = exp(tau)-1;
  const REAL psi2 = exp(pow(tau,2))-1;
  REAL k_a, k_aa, k_r, k_rr, k_i, k_ii, pid, pr, Dpr, pa;
  pid = 0.01*T*rho*Rf; //[bar] Ideal contauibution to pressure
  Dpr = 0.01*T*rho*Rf*(delta*F.ar_d(tau, delta) - delta*tau*F.ar_dt(tau, delta)); //[bar] Residual repulsive pressure
  pr = pid + Dpr; //[bar] Repulsive pressure
  pa = 0.01*F.P(T, rho) - pr; //[bar] Attauactive pressure
  k_a = (F.a_0 + F.a_1*psi1 + F.a_2*psi2)*tau;
  k_aa = (F.A_0 + F.A_1*psi1 + F.A_2*psi2)*pow(tau,3);
  k_r = (F.b_0 + F.b_1*psi1 + F.b_2*psi2)*tau;
  k_rr = (F.B_0 + F.B_1*psi1 + F.B_2*psi2)*pow(tau,3);
  k_i = (F.c_0 + F.c_1*psi1 + F.c_2*tau)*tau;
  k_ii = (F.C_0 + F.C_1*psi1 + F.C_2*psi2)*pow(tau,3);
  return F.eta0(T) + (k_i*pid + k_r*Dpr + k_a*pa + k_ii*pow(pid,2) + k_rr*pow(Dpr,2) + k_aa*pow(pa,2))*1e-3;
}

REAL cp0(const REAL& T)
{
  const REAL tau = F.Tc/T;
  const REAL Rf = R/F.MW;
  return Rf*(-0.629789+7.292937*pow(tau,-0.5)+5.154411*pow(tau,-0.75));
}

REAL k(const REAL& T, const REAL& rho)
{
  const REAL tau = T/F.Tc;
  const REAL delta = rho/F.rhoc;
  REAL k_crit = 0;
  
  if((tau>0.85) || (delta>0.85))
  {
    REAL arc = 1 - F.ac_10*acosh(1+F.ac_11*pow(pow((1-tau),2),F.ac_12));
    k_crit = (delta*exp((-pow(delta,F.ac_1)/F.ac_1)-pow((F.ac_2*(tau-1)),2)-pow(F.ac_3*(delta-1),2)))/ \
             pow((pow(pow((1-(1.0/tau)+F.ac_4*pow(pow((delta-1),2),(1/(2*F.ac_5)))),2),F.ac_6) + \
             pow(F.ac_7*pow(delta-arc,2),F.ac_8)),F.ac_9);
  }

  return F.k_c*(F.n_1*pow(tau,F.g_1)*pow(delta,F.h_1) + \
                    F.n_2*pow(tau,F.g_2)*pow(delta,F.h_2) + \
                    F.n_3*pow(tau,F.g_3)*pow(delta,F.h_3) + \
                    F.n_4*pow(tau,F.g_4)*pow(delta,F.h_4) + \
                    F.n_5*pow(tau,F.g_5)*pow(delta,F.h_5) + \
                    F.n_6*pow(tau,F.g_6)*pow(delta,F.h_6) + \
                    F.n_7*pow(tau,F.g_7)*pow(delta,F.h_7) + \
                    F.n_8*exp(-5*pow(delta,2))*pow(tau,g_8)*pow(delta,h_8) + \
                    F.n_c*k_crit);
}

} /* openprop */
#endif /* end of include guard: BASE_H */