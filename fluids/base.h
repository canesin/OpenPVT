//  OpenPROP - Open source properties library
//	Fabio Cesar Canesin <fabio.canesin@gmail.com>
//	MIT 2013 license
//	base.h - Implement default fluid struct for openprop.
#ifndef BASE_H
#define BASE_H

#include "../utils/utils.h"

namespace openprop
{
//Force definition of fluid as R134a - add basic testing in upper-level to verify
//if fluid calculated is the same of the fluid asking for calculations

	struct fluid
	{
		//Thermophysical Constants for R134a - default constants
		const REAL MW; //[g/mol] Molecular weight
		const REAL h0; //[kJ/kg] Reference enthalpy
		const REAL s0; //[kJ/kg*K] Reference enTopy at 273.15K
		const REAL Tb; //[K] Normal boiling point (at 1 atm)
		const REAL Tt; //[K] Triple point temperature
		const REAL Tc; //[K] Critical temperature
		const REAL Pc; //[MPa] Critical Pessure
		const REAL rhoc; //[kg/m3] Critical density
		const REAL omega; //Acentric factor
		const REAL sigma; //[nm] Minimum Leonard-Jones potential
		const REAL mu; //10^-30 [C*m] Dipole moment

		//Limit of validity for the fluid EoS
		const REAL T_max; //[K]
		const REAL T_min; //[K]
		const REAL P_min; //[kPa]
		const REAL P_max; //[kPa]

		// EoS constants for fluid
		const REAL SE1;
		const REAL SE2;
		const REAL SE3;
		const REAL SE4;
		const REAL SE5;
		const REAL SE6;
		const REAL SE7;
		const REAL SE8;
		const REAL SE9;
		const REAL SE10;
		const REAL SE11;
		const REAL SE12;
		const REAL SE13;
		const REAL SE14;

		//Friction theory constants
		const REAL a_0;
		const REAL a_1;
		const REAL a_2;
		const REAL b_0;
		const REAL b_1;
		const REAL b_2;
		const REAL c_0;
		const REAL c_1;
		const REAL c_2;
		const REAL d_0;
		const REAL d_1;
		const REAL d_2;
		const REAL d_3;
		const REAL A_0;
		const REAL A_1;
		const REAL A_2;
		const REAL B_0;
		const REAL B_1;
		const REAL B_2;
		const REAL C_0;
		const REAL C_1;
		const REAL C_2;

		//Thermal conductivity constants
		//Fluid Phase Equilibria 245 (2006) 37–51
		//G. Scalabrin et al.
		const REAL k_c;
		const REAL g_1;
		const REAL g_2;
		const REAL g_3;
		const REAL g_4;
		const REAL g_5;
		const REAL g_6;
		const REAL g_7;
		const REAL g_8;
		const REAL h_1;
		const REAL h_2;
		const REAL h_3;
		const REAL h_4;
		const REAL h_5;
		const REAL h_6;
		const REAL h_7;
		const REAL h_8;
		const REAL n_1;
		const REAL n_2;
		const REAL n_3;
		const REAL n_4;
		const REAL n_5;
		const REAL n_6;
		const REAL n_7;
		const REAL n_8;
		const REAL n_c;
		const REAL ac_1;
		const REAL ac_2;
		const REAL ac_3;
		const REAL ac_4;
		const REAL ac_5;
		const REAL ac_6;
		const REAL ac_7;
		const REAL ac_8;
		const REAL ac_9;
		const REAL ac_10;
		const REAL ac_11;
		const REAL ac_12;
	} /* fluid */
} /* openprop */
#endif /* end of include guard: BASE_H */