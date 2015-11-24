#ifndef _BIVARIATE_POLYNOMIAL_
#define  _BIVARIATE_POLYNOMIAL_

class BivariatePolynomial
{
private:
    int ** MatrixRepresentation;
    int matrixDimension;
    int degy,degx;
    void Parse(char str[]);
    void InitDegs();
public:
    int ** getMatrix()
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
    BivariatePolynomial(int systemDegree);
    void Print();                                               //print contents
    ~BivariatePolynomial();
    int getCoefFromPolynomial(int & pos, char * polstr);
    void ParsePolynomialBody(int & i, char * polstr, char first, char second,int & firstCoord, int & secondCoord, int & choice , int coef);
    double backSubstitute(double x, double y);

    double exp(double num, double power);

};

#endif