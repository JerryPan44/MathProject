#include "../../SylvesterMatrix.h"
#include "../../Polynomial.h"
#include "../../BivariatePolynomial.h"
#include "../../ChangeOfVariableCoefficients.h"
#include <gtest/gtest.h>

class SylvesterMatrixTests : public ::testing::Test
{
protected:
    SylvesterMatrix * SM;
    BivariatePolynomial * Bp1, * Bp2;
    virtual void SetUp() {
        char p1[20];
        strcpy(p1, "4+3*x^2+2*x*y");
        char p2[20];
        strcpy(p2, "3+4*x^2+5*x*y");
        Bp1 = new BivariatePolynomial(p1, 2);
        Bp1->Print();
        Bp2 = new BivariatePolynomial(p2, 2);
        Bp2->Print();
        SM = new SylvesterMatrix(Bp1, Bp2);
    }


};
TEST_F(SylvesterMatrixTests, testInitPolynomialPower)
{
    int p[2];
    Polynomial ** powerUp1;         //matrix with all the powers of t1*z+t2
    p[0] = 2;
    p[1] = 3;
    SM->initPowerUps(powerUp1, p, 4, 2);
    int * powerTest = powerUp1[3]->getPolynomial();
    EXPECT_EQ(8, powerTest[0]);
    EXPECT_EQ(36, powerTest[1]);
    EXPECT_EQ(54, powerTest[2]);
    EXPECT_EQ(27, powerTest[3]);
    p[0] = 4;
    p[1] = 8;
    Polynomial ** powerUp2;
    SM->initPowerUps(powerUp2, p, 4, 2);
    powerTest = powerUp2[3]->getPolynomial();
    EXPECT_EQ(64, powerTest[0]);
    EXPECT_EQ(384, powerTest[1]);
    EXPECT_EQ(768, powerTest[2]);
    EXPECT_EQ(512, powerTest[3]);
}
TEST_F(SylvesterMatrixTests, TestChangeOfVariableBody)
{
//    SM->Print();
    SM->coefs = new ChangeOfVariableCoefficients(1, 2, 3, 4);

    EXPECT_EQ(true, SM->changeOfVariableBody());
//    SM->Print();

}