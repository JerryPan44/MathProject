#ifndef _BIVARIATE_POLYNOMIAL_
#define  _BIVARIATE_POLYNOMIAL_
#include "Polynomial.h"
class BivariatePolynomial
{
private:
    double ** MatrixRepresentation;
    int matrixDimension;
    int degy,degx;
    void Parse(char str[]);
    void InitDegs();
public:
    double ** getMatrix()
    {
        return this->MatrixRepresentation;
    }
    int getMatrixDimension()
    {
        return this->matrixDimension;
    }

    int getdegy()
    {
        return degy;
    }
    int getdegx()
    {
        return degx;
    }
    bool consume(char * , int & , int );
    BivariatePolynomial(char polstr[], int systemDegree);
    BivariatePolynomial(int sD, double ** P);
    BivariatePolynomial(int systemDegree);
    void Print();                                               //print contents
    ~BivariatePolynomial();
    int getCoefFromPolynomial(int & pos, char * polstr);
    void ParsePolynomialBody(int & i, char * polstr, char first, char second,int & firstCoord, int & secondCoord, int & choice , double coef);
    double backSubstituteXandY(double x, double y);
    Polynomial * backSubstitute(double , char);
    double exp(double num, int power);
    double getDoubleCoefFromPolynomial(int & pos, char * polstr);
    int getIntCoefFromPolynomial(int & pos, char * polstr);
};

#endif
