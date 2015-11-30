#include <iostream>
#include "SylvesterMatrix.h"
#include "BivariatePolynomial.h"
#include <cstdlib>
#include <time.h>
using namespace std;

int SylvesterMatrix::max(int d1, int d2)							//Find max
{
    if (d1>d2)
        return d1;
    else
        return d2;
}

bool SylvesterMatrix::changeOfVariable()
{								//Pick 4 random t(i)
    int t1, t2, t3, t4;
    getRandom(t1);
    getRandom(t2);
    getRandom(t3);
    getRandom(t4);
    this->coefs = new ChangeOfVariableCoefficients(t1, t2, t3, t4);
    this->changeOfVariableBody();
}

bool SylvesterMatrix::changeOfVariableBody()
{
    int p[2];
    p[0] = coefs->t1;
    p[1] = coefs->t2;
    Polynomial ** powerUp1;									//matrix with all the powers of t1*z+t2
    initPowerUps(powerUp1, p, this->hiddenDeg + 1, 2);
    p[0] = coefs->t3;
    p[1] = coefs->t4;
    								//create list for all the powerups and pass this as an argument for the changeofVariable
    Polynomial ** powerUp2;         								//matrix with all the powers of t3*z+t4
    initPowerUps(powerUp2, p, this->hiddenDeg + 1, 2);
    for (int i = 0; i < this->RowDimension; ++i) {
        for (int j = 0; j < this->ColDimension; ++j) {
            this->matrix[i][j]->changeOfVariable(powerUp1, powerUp2, this->hiddenDeg);
        }
    }
    cleanPowerUps(powerUp1, this->hiddenDeg + 1);
    cleanPowerUps(powerUp2, this->hiddenDeg + 1);
    return true;
}
void SylvesterMatrix::cleanPowerUps(Polynomial **& p, int size)
{
    for (int i = 0; i < size; ++i) {
        delete p[i];
    }
    delete []p;
}
bool SylvesterMatrix::initPowerUps(Polynomial **&p, int * m, int sizeOfP, int sizeOfM)
{
    p = new Polynomial*[sizeOfP];
    for (int i = 0; i < sizeOfP; ++i) {
        p[i] = new Polynomial(sizeOfM - 1, m);
        p[i]->powerUp(i);
    }
    return true;

}
void SylvesterMatrix::getRandom(int & t)							//Random number
{
    srand (time(NULL)+rand()%10000);
    t = (rand()%100) - 50;
}
SylvesterMatrix::SylvesterMatrix(BivariatePolynomial * Bp1, BivariatePolynomial * Bp2)
{
    int maxDegx = max(Bp1->getdegx(),Bp2->getdegx());
    int maxDegy = max(Bp1->getdegy(),Bp2->getdegy());
    if(maxDegx>=maxDegy) {									//if maxdegx > maxdegy hide y
        d0 = Bp1->getdegx();
        d1 = Bp2->getdegx();
        cout<<"Hidden Variable is Y!"<<endl;
        this->hiddenVariable = 'y';
        this->hiddenDeg = maxDegy;
        this->initMatrixWithHiddenY(Bp1,Bp2);
    }
    else {											//else hide x
        d0 = Bp1->getdegy();
        d1 = Bp2->getdegy();
        this->hiddenDeg = maxDegx;
        this->hiddenVariable = 'x';
        cout<<"Hidden Variable is X!"<<endl;
        this->initMatrixWithHiddenX(Bp1,Bp2);
    }
}
void SylvesterMatrix::Print()									//Print Syl. Matrix
{
    cout<<"******"<<endl;
    cout<<"Sylvester matrix dimensions : "<<this->RowDimension<<"x"<<this->ColDimension<<endl;
    for (int i = 0; i < RowDimension; ++i) {
        for (int j = ColDimension-1; j >=0; --j) {
            this->matrix[i][j]->Print();
        }
        cout<<endl;
    }
    cout<<"******"<<endl;
}
void SylvesterMatrix::allocMatrix()
{
    RowDimension = d0 + d1;
    ColDimension = RowDimension;
    this->matrix = new Polynomial ** [RowDimension];
    for (int i = 0; i < ColDimension; ++i) {
        this->matrix[i] = new Polynomial * [ColDimension];
    }

    int * temp = new int[1];
    temp[0] = 0;
    for (int i = 0; i < RowDimension; ++i) {
        for (int j = 0; j < ColDimension; ++j) {
            this->matrix[i][j] = new Polynomial(0, temp);
        }
    }
    delete []temp;

}
void SylvesterMatrix::assignXCol(BivariatePolynomial * Bp, int col, int * & resultM)		//return the column col in resultM of the BivariatePolynomial
{
    int degy = Bp->getdegy();
    int ** matrix = Bp->getMatrix();
    resultM = new int[degy+1];
    int k=0;
    for (int i = 0; i < degy + 1; ++i) {
        resultM[k] = matrix[i][col];
        ++k;
    }
    //cout<<endl;
}

void SylvesterMatrix::assignYRow(BivariatePolynomial * Bp, int row, int * & resultM)		//return the row row in resultM of the BivariatePolynomial
{
    int degx = Bp->getdegx();
    int ** matrix = Bp->getMatrix();
    resultM = new int[degx+1];
    int k=0;
    for (int i = 0; i < degx + 1; ++i) {
        resultM[k] = matrix[row][i];
        ++k;
    }
    //cout<<endl;
}


void SylvesterMatrix::initMatrixWithHiddenY(BivariatePolynomial *Bp1, BivariatePolynomial *Bp2)         //Make the sylvester matrix with x as the hidden variable
{
    allocMatrix();
    int * temp, *ZeroTemp;
    int degy2=Bp2->getdegy(),degy1=Bp1->getdegy(),degx1 = Bp1->getdegx(), degx2 = Bp2->getdegx();
    int totalDeg = Bp1->getMatrixDimension();

    assignZero(degy1, ZeroTemp);
    for (int i = d1 - 1; i >= 0; --i) {
        int until = (d1-1)-i;
        for (int j = 0; j < ColDimension - until; ++j) {						//if j column beyond the degx assign zero
            if(j>degx1)
            {
                this->matrix[i][j + until] = new Polynomial(degy1, ZeroTemp);
                continue;
            }
            assignXCol(Bp1, j, temp);
            this->matrix[i][j + until] = new Polynomial(degy1, temp);
            delete []temp;
        }
    }
    delete []ZeroTemp;

    assignZero(degy2, ZeroTemp);
    for (int i = d1 + d0 - 1; i >= d1; --i) {
        int until = (d1+d0-1)-i;
        for (int j = 0; j < ColDimension - until; ++j) {						//if j column beyond the degx assign zero
            if(j>degx2)
            {
                this->matrix[i][j + until] = new Polynomial(degy2, ZeroTemp);
                continue;
            }
            assignXCol(Bp2, j, temp);
            this->matrix[i][j + until] = new Polynomial(degy2, temp);
            delete []temp;
        }
    }
    delete []ZeroTemp;


}
void SylvesterMatrix::assignZero(int deg, int*& temp)
{
    temp = new int[deg+1];
    for (int i = 0; i < deg + 1; ++i) {
        temp[i] = 0;
    }
}

void SylvesterMatrix::initMatrixWithHiddenX(BivariatePolynomial *Bp1, BivariatePolynomial *Bp2)		//Make the sylvester matrix with x as the hidden variable
{
    allocMatrix();
    int * temp, *ZeroTemp;
    int degy2=Bp2->getdegy(),degy1=Bp1->getdegy(),degx1 = Bp1->getdegx(), degx2 = Bp2->getdegx();
    int totalDeg = Bp1->getMatrixDimension();

    assignZero(degy1, ZeroTemp);
    for (int i = d1 - 1; i >= 0; --i) {
        int until = (d1-1)-i;
        for (int j = 0; j < ColDimension - until; ++j) {						//if j column beyond the degx assign zero
            if(j>degy1)
            {
                this->matrix[i][j + until] = new Polynomial(degx1, ZeroTemp);
                continue;
            }
            assignYRow(Bp1, j, temp);
            this->matrix[i][j + until] = new Polynomial(degx1, temp);
            delete []temp;
        }
    }
    delete []ZeroTemp;

    assignZero(degy2, ZeroTemp);
    for (int i = d1 + d0 - 1; i >= d1; --i) {
        int until = (d1+d0-1)-i;
        for (int j = 0; j < ColDimension - until; ++j) {						//if j column beyond the degx assign zero
            if(j>degy2)
            {
                this->matrix[i][j + until] = new Polynomial(degx2, ZeroTemp);
                continue;
            }
            assignYRow(Bp2, j, temp);
            this->matrix[i][j + until] = new Polynomial(degx2, temp);
            delete []temp;
        }
    }
    delete []ZeroTemp;
}

void SylvesterMatrix::multiply(MyMatrix *m, SylvesterMatrix*& result)			//multiply the sylvester matrix with m and store the result in result
{
    if(m == NULL)
        return;
    if(matrix == NULL)
        return;

    double ** matrix2 = m->getMatrix();
    Polynomial *** p = new Polynomial**[this->RowDimension];				//the result has the same number of rows as the sylvester matrix
    for (int i = 0; i < this->RowDimension; ++i) {
        p[i] = new Polynomial*[m->getColDimension()];					//and the same number of columns as m
        int rowMaxDegree = this->getMatrixRowMaxDegree(i);
        for (int j = 0; j < m->getColDimension(); ++j) {
            p[i][j] = new Polynomial(rowMaxDegree);		//the polynomials degree is going to be the same (of <=) as the max degree of the row of the sylvester
        }
    }
    Polynomial * temp;

    for (int k = 0; k < m->getColDimension(); ++k) {
        for (int i = 0; i < this->RowDimension; ++i) {
            for (int j = 0; j < this->ColDimension; ++j) {
                this->matrix[i][j]->multiply(temp, matrix2[j][k]);			//actual multiplication
                p[i][k]->Add(temp);					//add the result of the current column * row to the corresponding position of the final matrix
                delete temp;
            }
        }
    }
    result = new SylvesterMatrix(p, m->getColDimension(), this->RowDimension);
}
int SylvesterMatrix::getMatrixRowMaxDegree(int Row)					//get the maximum polynomial degree of the given row
{
    int max = 0;
    for (int i = 0; i < this->ColDimension; ++i) {
        if(this->matrix[Row][i]->getDegree()>max)
            max = this->matrix[Row][i]->getDegree();
    }
    return max;
}
SylvesterMatrix::SylvesterMatrix(Polynomial *** p, unsigned int cD, unsigned int rD)
{
    this->matrix = p;
    this->ColDimension = cD;
    this->RowDimension = rD;
}
SylvesterMatrix::SylvesterMatrix(SylvesterMatrix & SM)
{
    this->ColDimension = SM.getColDimension();
    this->RowDimension = SM.getRowDimension();
    this->hiddenDeg = SM.getHiddenDeg();
    this->matrix = new Polynomial**[this->RowDimension];
    for (int i = 0; i < this->RowDimension; ++i) {
        this->matrix[i] = new Polynomial*[this->ColDimension];
    }
    Polynomial *** p = SM.getMatrix();
    for (int j = 0; j < this->RowDimension; ++j) {
        for (int i = 0; i < this->ColDimension; ++i) {
            this->matrix[i][j] = new Polynomial(p[i][j]->getDegree(), p[i][j]->getPolynomial());
        }
    }

}
SylvesterMatrix::~SylvesterMatrix()
{
    for (int i = 0; i < RowDimension; ++i) {
        for (int j = 0; j < ColDimension; ++j) {
            delete this->matrix[i][j];
        }
        delete []this->matrix[i];
    }
    delete []this->matrix;
}
