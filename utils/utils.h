//  OpenPROP - Open source properties library
//	Fabio Cesar Canesin <fabio.canesin@gmail.com>
//	MIT 2013 license
//	utils.h - Utitilities for OpenPROP

#ifndef UTILS_H
#define UTILS_H

#define REALSMALL 1E-14
#define REALBIG 1E+14

//REAL is the definition for the precision in the algorithms of the library
//this forces the REAL, utils should prevail over local definitions
#define REAL double

//Define Gas constant - J mol-1 K-1
//NIST 2010 CODATA recommended value of R
#ifndef R
#define R 8.3144621
#endif

#ifndef PI
#define PI 3.141592653589793
#endif

#ifndef e
#define e 2.718281828459045
#endif

//Central diff function, makes two function calls, O(h^2)
REAL diff(const REAL h, const REAL x, REAL (*func)(const REAL))
{
	// diff = f(x + h) - f(x -h)/2h + O(h^2)
	return ((*func)(x + h) - (*func)(x - h))/(2.0*h + REALSMALL);
}

//Forward diff function, receive funcx making only one function call, O(h)
REAL diff_forward(const REAL h, const REAL x, const REAL funcx, REAL (*func)(const REAL))
{
	// diff = f(x + h) - f(x)/h + O(h)
	return ((*func)(x + h) - funcx)/(h + REALSMALL);
}

//Backward diff function, receive funcx making only one function call, O(h)
REAL diff_backward(const REAL h, const REAL x, const REAL funcx, REAL (*func)(const REAL))
{
	// diff = f(x) - f(x - h)/h + O(h)
	return (funcx - (*func)(x - h))/(h + REALSMALL);
}

//Central second diff function, receive funcx making only two function calls, O(h**2)
REAL diff2(const REAL h, const REAL x, const REAL funcx, REAL (*func)(const REAL))
{
	// diff2 = f(x +h) +2*f(x) - f(x - h)/h**2 + O(h**2)
	return ((*func)(x + h) + 2*funcx - (*func)(x - h))/(pow(h,2)s + REALSMALL);
}

//	Possitive unidimensional Newton-Raphson-Bissection algorithm
//	Use numerical forward differencing for speed (reduces number of function calls)
//	Use a initial guess when calling PNRB, like ideal gas equation solution for variable
REAL PNRB(REAL guess, REAL (*func)(const REAL))
{
	REAL fx = (*func)(guess);
	//Early-abort, check if guess is already root for equation
	if (fx < 1E-08 && guess > 0.0) return guess;
	//Lower limit for real space. Searching for a non-negative root
	REAL a = REALSMALL;
	//Upper limit for real space in REAL type
	REAL b = REALBIG;
	REAL fa = (*func)(a);
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
		guess = guess - fx/diff_forward(1.0E-8, guess, fx, func);
		//If the result is outside braket use bissection step (for stability)
		if ((b-guess)*(guess-a) < 0.0){
			guess = guess + 0.5*(b-a);
		}
		//Compute new value for function
		fx = (*func)(guess);
		//Check if guess is root for equation
		if (fx < 1E-08 && guess > 0.0){
			return guess;
		} 
		//Increase interation count
		count +=1;
	//Limit the iterations to 101, should converge much faster (like 3~5 iters)	
	} while(count < 101);
	return guess;
}

#endif /* end of include guard: UTILS_H */
