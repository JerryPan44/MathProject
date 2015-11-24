#ifndef _SYLVESTER_MATRIX_
#define _SYLVESTER_MATRIX_
#include "BivariatePolynomial.h"
#include "Polynomial.h"
#include "Matrix.h"
#include "ChangeOfVariableCoefficients.h"


class SylvesterMatrix
{
    unsigned int RowDimension, ColDimension;
    unsigned int d0,d1;
    int hiddenDeg;
    Polynomial *** matrix;
    char hiddenVariable;
    void assignXCol(BivariatePolynomial * Bp, int power, int * & resultM);
    void assignYRow(BivariatePolynomial * Bp, int row, int * & resultM);
    void initMatrixWithHiddenX(BivariatePolynomial * Bp1, BivariatePolynomial * Bp2);
    void initMatrixWithHiddenY(BivariatePolynomial * Bp1, BivariatePolynomial * Bp2);
    void cleanPowerUps(Polynomial **&, int);
public:
    ChangeOfVariableCoefficients * coefs;
    static bool initPowerUps(Polynomial **& ,  int * , int, int);
    ChangeOfVariableCoefficients * getCoefs()
    {
        return this->coefs;
    }
    int getHiddenDeg()
    {
        return hiddenDeg;
    }
    unsigned int getRowDimension()
    {
        return  RowDimension;
    }
    unsigned int getColDimension()
    {
        return  ColDimension;
    }
    char getHiddenVariable()
    {
        return  hiddenVariable;
    }
    Polynomial *** getMatrix()
    {
        return this->matrix;
    }
    int max(int, int);
    SylvesterMatrix(BivariatePolynomial * Bp1, BivariatePolynomial * Bp2);
    SylvesterMatrix(Polynomial ***, unsigned int, unsigned int);
    SylvesterMatrix(SylvesterMatrix & SM);
    bool changeOfVariableBody();
    void allocMatrix();
    bool changeOfVariable();
    void Print();
    void multiply(MyMatrix * m, SylvesterMatrix*& result);
    int getMatrixRowMaxDegree(int Row);
    void assignZero(int deg, int*& temp);
    void getRandom(int &);
    ~SylvesterMatrix();
};

#endif
