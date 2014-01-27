//  OpenPROP - Open source properties library
//	Fabio Cesar Canesin <fabio.canesin@gmail.com>
//	MIT 2013 license
//	utils.h - Utitilities for OpenPROP

#ifndef UTILS_H
#define UTILS_H

#define REALSMALL 1E-14
#define REALBIG 1E+14
#define TOLDEFAULT 1.5E-07

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

//Helper function to see if value is in tolerance
int is_equal(double a, double b, const double epsilon = TOLDEFAULT);
//Central diff function, makes two function calls, O(h^2)
REAL diff(const REAL x, REAL (*func)(const REAL), const REAL h = TOLDEFAULT);
//Forward diff function, receive funcx making only one function call, O(h)
REAL diff_forward(const REAL x, const REAL funcx, REAL (*func)(const REAL), const REAL h = TOLDEFAULT);
//Backward diff function, receive funcx making only one function call, O(h)
REAL diff_backward(const REAL x, const REAL funcx, REAL (*func)(const REAL), const REAL h = TOLDEFAULT);
//Central second diff function, receive funcx making only two function calls, O(h**2)
REAL diff2(const REAL x, const REAL funcx, REAL (*func)(const REAL), const REAL h = TOLDEFAULT);
//	Possitive unidimensional Newton-Raphson-Bissection algorithm
//	Use numerical forward differencing for speed (reduces number of function calls)
//	Use a initial guess when calling PNRB, like ideal gas equation solution for variable
REAL PNRB(REAL guess, REAL (*func)(const REAL), const REAL h = 0.1*TOLDEFAULT);

#endif /* end of include guard: UTILS_H */
