#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "../utils/utils.h"

REAL x1(const REAL x);
REAL x3x2(const REAL x);
REAL x1_0001(const REAL x);
REAL x2_0025(const REAL x);

int main ()
{
	REAL res;
	REAL tol = 1E-5;
	//Testing central diff
	printf("Central diff\n");
	res = diff(1, *x1);
	assert(is_equal(1, res, tol));
	printf("%.8f \n", res);
	res = diff(1, *x3x2);
	assert(is_equal(5, res, tol));
	printf("%.8f \n", res);
	
	//Testing forward diff
	printf("Forward diff\n");
	res = diff_forward(1, x1(1), *x1);
	assert(is_equal(1, res, tol));
	printf("%.8f \n", res);
	res = diff_forward(1, x3x2(1), *x3x2);
	assert(is_equal(5, res, tol));
	printf("%.8f \n", res);

	//Testing backward diff
	printf("Backward diff\n");
	res = diff_backward(1, 1, *x1);
	assert(is_equal(1, res, tol));
	printf("%.8f \n", res);
	res = diff_backward(1, 2, *x3x2);
	assert(is_equal(5, res, tol));
	printf("%.8f \n", res);

	//Testing central diff2
	printf("Central diff2\n");
	res = diff2(1, 1, *x1);
	assert(is_equal(0, res, tol));
	printf("%.8f \n", res);
	res = diff2(1, 2, *x3x2);
	assert(is_equal(8.0, res, 0.01));
	printf("%.8f \n", res);

	//Testing PNRB
	printf("Root by PNRB\n");
	res = PNRB(1513123.131533, *x1_0001);
	assert(is_equal(1E-4, res, tol));
	printf("%.8f \n", res);
	res = PNRB(1514563.123124, *x2_0025);
	assert(is_equal(0.0005, res, tol));
	printf("%.8f \n", res);

	return 0;
}

REAL x1_0001(const REAL x)
{
	return x -1E-4;
}

REAL x2_0025(const REAL x)
{
	return x*x -0.00000025;
}

REAL x1(const REAL x)
{
	return x;
}

REAL x3x2(const REAL x)
{
	return x*x*x + x*x;
}