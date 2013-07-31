//  OpenPROP - Open source properties library
//	Fabio Cesar Canesin <fabio.canesin@gmail.com>
//	MIT 2013 license
//	R134a.h - R134a (1,1,1,2-tetrafluoroethane) fluid file.

#include "../utils/utils.h"

//Force definition of fluid as R134a - add basic testing in upper-level to verify
//if fluid calculated is the same of the fluid asking for calculations
#define fluid R134a

struct fluid{
	REAL Tc = 374.083; //Kelvin - critical temperature
	REAL pc = 4.048; //MPa - critical pressure
	REAL rhoc = 509; //kgm-3 - critical density
	REAL Tt = 169.85; //Kelvin - triple temperature
	REAL M = 0.102031; //kg mol-1 - molecular mass
}

//Reference for: p_sat, rhol_sat, rhog_sat and fluid properties
//A fundamental equation of state for 1,1,1,2-tetrafluoroethane with an intermolecular
//potential energy background and reliable ideal-gas properties
//I Made Astina, Haruki Sato
//Fluid Phase Equilibria 221 (2004) 103–111

//Saturation pressure
REAL p_sat(REAL T){
	REAL t = T/fluid.Tc; //Reduced temperature
	return fluid.pc * exp((1/t)*(
		-7.628631*pow(1-t,1) + \
		1.734314*pow(1-t,1.5) + \
		−2.591127*pow(1-t,2.5) + \
		−3.290324*pow(1-t,5)
		));
}

//Saturated gas density
REAL rhog_sat(REAL T){
	REAL t = T/fluid.Tc; //Reduced temperature
	return fluid.pc * exp((1/t)*(
		-7.628631*pow(1-t,1) + \
		1.734314*pow(1-t,1.5) + \
		−2.591127*pow(1-t,2.5) + \
		−3.290324*pow(1-t,5)
		));
}

//Saturated liquid density
REAL rhol_sat(REAL T){
	REAL t = T/fluid.Tc; //Reduced temperature
	return fluid.pc * ((1/t)*(
		1.8689487*pow(1-t,1/3) + \
		0.6531694*pow(1-t,2/3) + \
		0.22950915*pow(1-t,1) + \
		0.7436063*pow(1-t,5)
		)+1);
}