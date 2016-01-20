//#include <limits.h>
#include "../../Polynomial.h"
#include "../../SylvesterMatrix.h"
#include <gtest/gtest.h>

class PolynomialTests : public ::testing::Test
{

protected:
    Polynomial * Pol, * Pol2;
    virtual void SetUp() {
        double p[3];
        p[0] = 5;
        p[1] = 3;
        p[2] = 2;
        Pol = new Polynomial(2, p);
        p[0] = 2;
        p[1] = 4;
        p[2] = 3;
        Pol2 = new Polynomial(2, p);
    }


};

TEST_F(PolynomialTests, testPowerUp1)
{
    Pol->powerUp(2);
    double * pol = Pol->getPolynomial();
    EXPECT_EQ(4, Pol->getDegree());
    EXPECT_EQ(4, pol[4]);
    EXPECT_EQ(12, pol[3]);
    EXPECT_EQ(29, pol[2]);
    EXPECT_EQ(30, pol[1]);
    EXPECT_EQ(25, pol[0]);

}

TEST_F(PolynomialTests, testPowerUp2)
{
    Pol->powerUp(3);
    double * pol = Pol->getPolynomial();
    EXPECT_EQ(6, Pol->getDegree());
    EXPECT_EQ(8, pol[6]);
}


TEST_F(PolynomialTests, testMultiplyPolynomials)
{
    Polynomial * result;
    Pol->multiplyPolynomials(Pol, Pol2, result);
    EXPECT_EQ(4, result->getDegree());
    double  * p = result->getPolynomial();
    EXPECT_EQ(6, p[4]);
    EXPECT_EQ(17, p[3]);
    EXPECT_EQ(31, p[2]);
    EXPECT_EQ(26, p[1]);
    EXPECT_EQ(10, p[0]);
}

//TEST_F(PolynomialTests, testMultiplyPolynomials)
//{
//Polynomial * result;
//Pol->multiplyPolynomials(Pol, Pol2, result);
//EXPECT_EQ(4, result->getDegree());
//int  * p = result->getPolynomial();
//EXPECT_EQ(6, p[4]);
//EXPECT_EQ(17, p[3]);
//EXPECT_EQ(31, p[2]);
//EXPECT_EQ(26, p[1]);
//EXPECT_EQ(10, p[0]);
//}

TEST_F(PolynomialTests, testChangeOfVariable)
{
    Polynomial * result;
    double testM[3][3];
    Polynomial ** p1, **p2;
    p1 = new Polynomial*[2];
    Polynomial ** powerUp1;         //matrix with all the powers of t1*z+t2
    double p[2];
    p[0] = 2;
    p[1] = 3;
    SylvesterMatrix::initPowerUps(powerUp1, p, 3, 2);
    p[0] = 3;
    p[1] = 4;
    Polynomial ** powerUp2;
    SylvesterMatrix::initPowerUps(powerUp2, p, 3, 2);
    Pol->changeOfVariable(powerUp1, powerUp2, 2);
    double * ret = Pol->getPolynomial();
    EXPECT_EQ(71, ret[0]);
    EXPECT_EQ(195, ret[1]);
    EXPECT_EQ(134, ret[2]);
}