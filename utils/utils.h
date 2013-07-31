#ifndef UTILS_H
#define UTILS_H

//  OpenPROP - Open source properties library
//	Fabio Cesar Canesin <fabio.canesin@gmail.com>
//	MIT 2013 license
//	utils.h - Utitilities for OpenPROP
#include <functional>
#include <limits>

//REAL is the definition for the precision in the algorithms of the library
//this forces the REAL, utils should prevail over local definitions
#define REAL double

//Define Gas constant - J mol-1 K-1
#ifndef R
#define R 8.314472
#endif

#ifndef PI
#define PI 3.141592653589793
#endif

#ifndef e
#define e 2.718281828459045
#endif


//Central difference function, makes two function calls, O(h^2)
REAL diff_central(const REAL& h, const REAL& x, const std::function<REAL (const REAL&)>& func)
{
	// diff = f(x + h) - f(x -h)/2h + O(h^2)
	return (func(x + h) - func(x - h))/(2.0*h + std::numeric_limits<REAL>::min());
}

//Foward difference function, receive funcx making only one function call, O(h)
REAL diff_foward(const REAL& h, const REAL& x, const REAL& funcx, const std::function<REAL (const REAL&)>& func)
{
	// diff = f(x + h) - f(x)/h + O(h)
	return (func(x + h) - funcx)/(h + std::numeric_limits<REAL>::min());
}

//	Possitive unidimensional Newton-Raphson-Bissection algorithm
//	Use numerical foward differencing for speed (reduces number of function calls)
//	Use a initial guess when calling PNRB, like ideal gas equation solution for variable
REAL PNRB(REAL guess, const std::function<REAL (const REAL&)>& func)
{
	REAL fx = func(guess);
	//Early-abort, check if guess is already root for equation
	if (fx < 1E-12 && guess > 0.0) return guess;
	//Lower limit for real space. Searching for a non-negative root
	REAL a = std::numeric_limits<REAL>::min();
	//Upper limit for real space in REAL type
	REAL b = std::numeric_limits<REAL>::max();
	REAL fa = func(a);
    int count = 0;
	do {
		//If root is in [a, guess] reduces search space around it.
		if (fx*fa < 0.0){
			b = guess;
		}
		//Else advance to guess position
		else { 
			a = guess;
		}
		//Choose new guess using Newton approximation
		guess = guess - fx/diff_foward(1.0E-5, guess, fx, func);
		//If the result is outside braket use bissection step (for stability)
		if ((b-guess)*(guess-a) < 0.0) guess = guess + 0.5*(b-a);
		//Compute new value for function
		fx = func(guess);
		//Check if guess is root for equation
		if (fx < 1E-12 && guess > 0.0) return guess;
		//Increase interation count
		count +=1;
	//Limit the iterations to 101, should converge much faster (like 4~6 iters)	
	} while(count < 101);
	return guess;
}

#endif /* end of include guard: UTILS_H */