#ifndef _SYLVESTER_POLYNOMIAL_
#define _SYLVESTER_POLYNOMIAL_

#include "SylvesterMatrix.h"
#include "Matrix.h"
#include "eigen/Eigen/Dense"
#include "eigen/Eigen/Core"
#include "eigen/Eigen/SVD"
class SylvesterPolynomial
{
    unsigned int degree, matrixDimensions;
    MyMatrix ** polynomial;
public:
    SylvesterPolynomial(int deg, int dimensions);
    void AssignPolynomial(int, Polynomial *** );
    bool SMatrixToSPolynomial(SylvesterMatrix * SM);
    void Print();
    unsigned int getMatrixDimensions()
    {
	    return this->matrixDimensions;
    }
    unsigned int getDegree()
    {
        return this->degree;
    }
    double ** getMd(int d)
    {
        if(d < 0 || d > this->matrixDimensions)
            return 0;
        return polynomial[d]->getMatrix();
    }
};

#endif
