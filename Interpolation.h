#ifndef _INTERPOLATION_
#define _INTERPOLATION_

#include <iostream>
#include "eigen/Eigen/Dense"
#include "BivariatePolynomial.h"

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
	Interpolation(unsigned int, unsigned int,MatrixXd);
	void ComputeM();
	double entry(unsigned int, unsigned int, unsigned int);
	int Check_k();
	double powerOf(double, unsigned int);
	BivariatePolynomial* find();
};

#endif
