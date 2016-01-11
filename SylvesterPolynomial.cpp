#include <iostream>
#include "SylvesterPolynomial.h"
using namespace std;
SylvesterPolynomial::SylvesterPolynomial(int deg, int dimensions)						//Create Polynomial for Sylvester Matrix
{
    this->polynomial = new MyMatrix *[deg + 1];
    for (int i = 0; i < deg + 1; ++i) {
        this->polynomial[i] = new MyMatrix(dimensions);								//Every polynomial is a MyMatrix
    }
    this->matrixDimensions = dimensions;
    this->degree = deg;
}


bool SylvesterPolynomial::SMatrixToSPolynomial(SylvesterMatrix * SM)						//Syl. Matrix to Syl. Polynomial 
{
    this->hiddenVar = SM->getHiddenVariable();
    if(SM->getColDimension() != this->matrixDimensions && SM->getRowDimension() != this->matrixDimensions)	//Dimensions don't match
        return  false;
    Polynomial *** PolMatrix = SM->getMatrix();
    for (int i = 0; i < this->degree + 1; ++i) {
        AssignPolynomial(i, PolMatrix);
    }
}
void SylvesterPolynomial::AssignPolynomial(int currDeg, Polynomial *** PolMatrix)				//Sylvester polynomial to matrix representation
{
    double ** matrix = this->polynomial[currDeg]->getMatrix();
    for (int i = 0; i < matrixDimensions; ++i) {
        for (int j = 0; j < matrixDimensions; ++j) {
            if(currDeg > PolMatrix[i][j]->getDegree())								//If there is no variable over max degree
                matrix[i][j]=0;											//fill matrix with zeros
            else
            {
                matrix[i][j] = PolMatrix[i][j]->getPolynomial()[currDeg];
            }
        }
    }
}

void SylvesterPolynomial::Print()										//Print Syl. Polynomial
{
    cout<<"*******"<<endl;
    cout<<"Sylvester Polynomial Degree : "<<this->degree<<endl;
    for (int i = 0; i < this->degree + 1; ++i) {
        if(i==0)
            cout<<" + ";
        else
            cout<<" + y^"<<i<<"*";
        cout<<"matrix"<<i;
    }
    cout<<endl;
    for (int i = 0; i < this->degree + 1; ++i) {
        cout<<"matrix"<<i<<" : "<<endl;
        this->polynomial[i]->Print();
    }
    cout<<"*******"<<endl;
}
