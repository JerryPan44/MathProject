#ifndef _INTERPOLATION_
#define _INTERPOLATION_

#include <iostream>
#include "eigen/Eigen/Dense"
#include "BivariatePolynomial.h"
#include "inttypes.h"
using namespace Eigen;
class Interpolation
{
	int d;
	int k;
	int deg;
	MatrixXd A;
	MatrixXd M;
	MatrixXd ker;
	double** P;
	int rank;
 public:
	Interpolation(uint64_t, MatrixXd);
	void ComputeM();
	double entry(uint64_t, uint64_t, uint64_t);
	int Check_k();
	double powerOf(double, uint64_t);
	BivariatePolynomial* find();
};

#endif
