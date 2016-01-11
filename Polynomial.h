#ifndef _POLYNOMIAL_
#define  _POLYNOMIAL_

class Polynomial
{
    double * matrixPolynomial;
    int degree;
public:
    Polynomial();
    Polynomial(int deg, double * matrix);
    Polynomial(int deg);
    Polynomial(Polynomial&);
    ~Polynomial()
    {
        delete []this->matrixPolynomial;
    }
    double * getPolynomial()
    {
        return matrixPolynomial;
    }
    int getDegree()
    {
        return this->degree;
    }
    void Print();
    void multiply(int);
    void multiply(Polynomial*& RetPol, int factor);
    void Add(Polynomial* temp);
    bool changeOfVariable(Polynomial ** p1, Polynomial ** p2, int deg);
    bool powerUp(int power);
    bool multiplyPolynomial(Polynomial * p2);
    static void multiplyPolynomials(Polynomial * , Polynomial * , Polynomial * & );
    double * computeAndGetRoots(int &);
};

#endif
