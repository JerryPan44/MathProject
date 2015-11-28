#include <iostream>
#include <cstdlib>
#include "Matrix.h"
using namespace std;
MyMatrix::MyMatrix(int dim)
{
    this->Rowdimension = dim;
    this->Coldimension = dim;
    this->matrix = new double*[this->Rowdimension];
    for (int i = 0; i < this->Rowdimension; ++i) {
        this->matrix[i] = new double[this->Rowdimension];
    }
}

void MyMatrix::Print()
{
    for (int i = 0; i < this->Rowdimension ; ++i) {
        for (int j = this->Coldimension - 1; j >= 0; --j) {
            cout << this->matrix[i][j]<<" ";
        }
        cout<<endl;
    }
}
MyMatrix::MyMatrix(int rowDim, int colDim)
{
    this->Rowdimension = rowDim;
    this->Coldimension = colDim;
    this->matrix = new double*[this->Rowdimension];
    for (int i = 0; i < this->Rowdimension; ++i) {
        this->matrix[i] = new double[this->Coldimension];
    }
}

void MyMatrix::generate()             //random matrix generation
{
    srand (time(NULL));
    for (int i = 0; i < Rowdimension; ++i) {
        for (int j = 0; j < Coldimension; ++j) {
            this->matrix[i][j] = (rand()%100)-50;
        }
    }
}