#include <math.h>
#include <stdio.h>
#include "utils.h"

//Helper function to see if value is in tolerance
int is_equal(double a, double b, const double epsilon)
{
    double c = a - b;
    return fabs(c) <= epsilon;
}

//Central diff function, makes two function calls, O(h^2)
REAL diff(const REAL x, REAL (*func)(const REAL), const REAL h)
{
	// diff = f(x + h) - f(x -h)/2h + O(h^2)
	return ((*func)(x + h) - (*func)(x - h))/(2.0*h + REALSMALL);
}

//Forward diff function, receive funcx making only one function call, O(h)
REAL diff_forward(const REAL x, const REAL funcx, REAL (*func)(const REAL), const REAL h)
{
	// diff = f(x + h) - f(x)/h + O(h)
	return ((*func)(x + h) - funcx)/(h + REALSMALL);
}

//Backward diff function, receive funcx making only one function call, O(h)
REAL diff_backward(const REAL x, const REAL funcx, REAL (*func)(const REAL), const REAL h)
{
	// diff = f(x) - f(x - h)/h + O(h)
	return (funcx - (*func)(x - h))/(h + REALSMALL);
}

//Central second diff function, receive funcx making only two function calls, O(h**2)
REAL diff2(const REAL x, const REAL funcx, REAL (*func)(const REAL), const REAL h)
{
	// diff2 = f(x +h) -2*f(x) +f(x - h)/h**2 + O(h**2)
	return ((*func)(x + h) -2.0*funcx +(*func)(x - h))/(h*h +REALSMALL*REALSMALL);
}

//	Possitive unidimensional Newton-Raphson-Bissection algorithm
//	Use numerical forward differencing for speed (reduces number of function calls)
//	Use a initial guess when calling PNRB, like ideal gas equation solution for variable
REAL PNRB(REAL guess, REAL (*func)(const REAL), const REAL h)
{
	REAL fx = (*func)(guess);
	//Early-abort, check if guess is already root for equation
	if (fx < h && guess > 0.0) return guess;
	//Lower limit for real space. Searching for a non-negative root
	REAL a = REALSMALL;
	//Upper limit for real space in REAL type
	REAL b = REALBIG;
	REAL fa = (*func)(a);
	unsigned short int count = 0;
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
		guess = guess - fx/diff_forward(guess, fx, func, h);
		
		//If the result is outside braket use bissection step (for stability)
		if ((b-guess)*(guess-a) < 0.0){
			guess = guess + 0.5*(b-a);
		}
		//Compute new value for function
		fx = (*func)(guess);
		//Check if guess is root for equation
		if (fx < h && guess > 0.0){
			return guess;
		} 
		//Increase interation count
		count +=1;
	//Limit the iterations to 1001, should converge much faster (like 5~10 iters)	
	} while(count < 1001);
	return guess;
}