#ifndef _INTERPOLATION_
#define _INTERPOLATION_

#include <iostream>
#include <fstream>
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
	int getRank()
	{
		return rank;
	}
	int getK()
	{
		return k;
	}
	Interpolation(uint64_t, MatrixXd);
	~Interpolation();
	MatrixXd & ComputeM();
	double entry(uint64_t, uint64_t, uint64_t);
	int Check_k();
	double powerOf(double, uint64_t);
	BivariatePolynomial * find(std::ofstream &);
};

#endif
