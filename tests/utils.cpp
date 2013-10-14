#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "../utils/utils.h"

REAL x1(const REAL x)
{
	return x;
}

REAL x3x2(const REAL x)
{
	return x*x*x + x*x;
}

int is_equal(double a, double b, const double epsilon)
{
    double c = a - b;
    return fabs(c) < epsilon;
}

int main ()
{
	REAL res;
	REAL tol = 1e-5;
	REAL tol2 = pow(tol, 2);
	//Testing central diff
	printf("Central diff\n");
	res = diff(tol, 1, *x1);
//	assert(is_equal(1, res, tol2));
	printf("%.8f \n", res);
	res = diff(tol, 1, *x3x2);
//	assert(is_equal(5, res, tol2));
	printf("%.8f \n", res);
	
	//Testing forward diff
	printf("Forward diff\n");
	res = diff_forward(tol, 1, x1(1), *x1);
//	assert(is_equal(1, res, tol2));
	printf("%.8f \n", res);
	res = diff_forward(tol, 1, x3x2(1), *x3x2);
//	assert(is_equal(5, res, 0.0001));
	printf("%.8f \n", res);

	//Testing backward diff
	printf("Backward diff\n");
	res = diff_backward(tol, 1, 1, *x1);
//	assert(is_equal(1, res, tol2));
	printf("%.8f \n", res);
	res = diff_backward(tol, 1, 2, *x3x2);
//	assert(is_equal(5, res, 0.0001));
	printf("%.8f \n", res);

	//Testing central diff2
	printf("Central diff2\n");
	res = diff2(tol, 1, 1, *x1);
//	assert(is_equal(0.0, res, tol2));
	printf("%.8f \n", res);
	res = diff2(tol, 1, 2, *x3x2);
//	assert(is_equal(7, res, 0.0001));
	printf("%.8f \n", res);


//	REAL diff2(const REAL& h, const REAL& x, const REAL& funcx, REAL (*func)(const REAL&));

	return 0;
}
