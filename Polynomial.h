#ifndef _POLYNOMIAL_
#define  _POLYNOMIAL_

class Polynomial
{
    int * matrixPolynomial;
    int degree;
public:
    Polynomial();
    Polynomial(int deg, int matrix[]);
    Polynomial(int deg);
    Polynomial(Polynomial&);
    ~Polynomial()
    {
        delete []this->matrixPolynomial;
    }
    int * getPolynomial()
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

};

#endif
