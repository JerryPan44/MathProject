#include <iostream>
#include "Polynomial.h"
#include <cmath>
#include <assert.h>
#include "SylvesterPolynomial.h"
#include "eigen/Eigen/Dense"
#include "eigen/Eigen/Core"
#include "eigen/Eigen/SVD"
#include "eigen/Eigen/Eigenvalues"
#include "ErrorMargin.h"
using namespace std;
using namespace Eigen;

Polynomial::Polynomial()								//Polynomila class constructor (nothing given)
{
    this->degree=-1;
}
Polynomial::Polynomial(int deg, double * matrix)
{
    this->degree = deg;
    this->matrixPolynomial = new double[deg + 1];
    for (int i = 0; i < deg+1; ++i) {
        this->matrixPolynomial[i] = matrix[i];
    }
}

Polynomial::Polynomial(int deg)								//Create Polynomial from given degree, zero matrix
{
    this->degree = deg;
    this->matrixPolynomial = new double[deg + 1];
    for (int i = 0; i < deg+1; ++i) {
        this->matrixPolynomial[i] = 0;
    }

}
void Polynomial::Print()								//Print Polynomial
{
    if(degree == 0)
    {

        cout<<" ( "<<0<<" ) ";
        return;
    }
    int allZero = 0;
    cout<<" ( ";
    for (int i = 0; i < degree + 1; ++i) {						//For all the positions of the matrix
        if(this->matrixPolynomial[i] != 0) {
            if (i == 0)
                cout << this->matrixPolynomial[i];
            else
            {

                if(this->matrixPolynomial[i]>0) {					//print coeficients of y^1
                    if(i == 1) {
                        cout << "+" << this->matrixPolynomial[i] << "*y";
                        continue;
                    }
                    cout << "+" << this->matrixPolynomial[i] << "*y^" << i;
                }
                else {
                    if(i == 1)
                    {
                        cout << this->matrixPolynomial[i] << "*y";
                        continue;
                    }
                    cout << this->matrixPolynomial[i] << "*y^" << i;
                }
            }
        }
        else
        {
            allZero++;
        }

        }
    if(allZero == degree + 1)								//if all the coefs are zero print 0
        cout<<"0";
    cout<<" ) ";
}

void Polynomial::multiply(Polynomial*& RetPol, int factor)				//multiply the polynomial with factor and assign/return the result in RetPol
{
    RetPol = new Polynomial(this->degree, this->matrixPolynomial);
    RetPol->multiply(factor);
}
void Polynomial::multiply(int factor)							//multiply polynomial with a factor(integer)
{
    for (int i = 0; i < this->degree + 1; ++i) {
        this->matrixPolynomial[i] *= factor;
    }
}
void Polynomial::Add(Polynomial * temp)							//add temp to this polynomial
{
    if(this->degree < temp->getDegree())
        return;
    for (int i = 0; i < temp->getDegree() + 1; ++i) {
        this->matrixPolynomial[i] += temp->matrixPolynomial[i];
    }
}


bool Polynomial::changeOfVariable(Polynomial ** p1, Polynomial ** p2, int deg)
{
    if(this->degree == 0)
        return true;
    double * newMatrix = new double [deg + 1];
    for (int k = 0; k < deg + 1; ++k) {
        newMatrix[k] = 0;
    }
    Polynomial * tempResult;
    for (int i = 0; i < this->degree + 1; ++i) {
        this->multiplyPolynomials(p1[i], p2[deg - i], tempResult);
        double * tempMatrix = tempResult->getPolynomial();
        for (int j = 0; j < tempResult->getDegree() + 1; ++j) {
            newMatrix[j] += this->matrixPolynomial[i] * tempMatrix[j];
        }
        delete tempResult;
    }
    delete [] this->matrixPolynomial;
    this->matrixPolynomial = newMatrix;
}
Polynomial::Polynomial(Polynomial &p)
{
    this->degree = p.getDegree();
    this->matrixPolynomial = new double [this->degree + 1];
    double * oldPol = p.getPolynomial();
    for (int i = 0; i < this->degree + 1; ++i) {
        this->matrixPolynomial[i] = oldPol[i];
    }
}
bool Polynomial::powerUp(int power)
{
    if(power < 0)
        return false;
    if(power == 0)
    {
        this->matrixPolynomial[0] = 1;
        for (int i = 1; i < this->degree + 1; ++i) {
            this->matrixPolynomial[i] = 0;
        }
    }
    Polynomial * self = new Polynomial(*this);
    for (int i = 0; i < power - 1; ++i) {
        this->multiplyPolynomial(self);
    }
    return true;
}
bool Polynomial::multiplyPolynomial(Polynomial * p2)
{
    int size1 = this->degree + 1;
    int size2 = p2->getDegree() + 1;
    double * newPol = new double [size1 + size2 - 1];
    for (int k = 0; k < size1 + size2 - 1; ++k) {
        newPol[k] = 0;
    }
    double * p2Matrix = p2->getPolynomial();
    for (int i = 0; i < size1; ++i) {
        for (int j = 0; j < size2; ++j) {
            //multiply
            newPol[i+j] += this->matrixPolynomial[i] * p2Matrix[j];
        }
    }
    delete []this->matrixPolynomial;
    this->matrixPolynomial = newPol;
    this->degree = this->degree + p2->getDegree();
}

void Polynomial::multiplyPolynomials(Polynomial * p1, Polynomial * p2, Polynomial * & resultPol)
{										//multiply two polynomials
    int size1 = p1->getDegree() + 1;
    int size2 = p2->getDegree() + 1;
    double * newPol = new double [size1 + size2 - 1];
    for (int k = 0; k < size1 + size2 - 1; ++k) {
        newPol[k] = 0;
    }
    double * p1Matrix = p1->getPolynomial();
    double * p2Matrix = p2->getPolynomial();
    for (int i = 0; i < size1; ++i) {
        for (int j = 0; j < size2; ++j) {
            //multiply
            newPol[i+j] += p1Matrix[i] * p2Matrix[j];
        }
    }
    resultPol = new Polynomial(p1->getDegree() + p2->getDegree(), newPol);
}

double * Polynomial::computeAndGetRoots(int & numOfRoots)
{
	MatrixXd companion = MatrixXd::Zero(this->degree, this->degree);
	double temp[this->degree];
//    cout<<"***polynomial***"<<endl;
//    this->Print();
//    cout<<"******"<<endl;
    assert(this->matrixPolynomial[this->degree] != 0);
	for(int i =0; i < this->degree; i++)
	{
		temp[i] = this->matrixPolynomial[i]/this->matrixPolynomial[this->degree];
	}
	for(int i =0; i < this->degree - 1; i++)
	{
		companion(i, i+1) = 1;
	}
	for(int i =0; i < this->degree; i++)
	{
		companion(this->degree -1, i) = -temp[i];
	}
    EigenSolver<MatrixXd> es(companion);				//Give created Companion matrix to eigen lib
    MatrixXcd eivals = es.eigenvalues();				//Eigenvalues
    double * roots;
    numOfRoots = 0;
//    cout<<"*********"<<endl;
//    cout<<eivals<<endl;
//    cout<<"*********"<<endl;
    for(int i =0; i < eivals.rows(); i++)
	{
		if(abs(eivals(i).imag()) < ERROR_MARGIN)
		{
			numOfRoots++;
		}
	}
	roots = new double[numOfRoots];
	for(int i =0; i < eivals.rows(); i++)
	{
		if(abs(eivals(i).imag()) < ERROR_MARGIN)
		{
			roots[i] = eivals(i).real();
		}
	}
	return roots;
}
