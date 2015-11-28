#include "../../BivariatePolynomial.h"
#include <gtest/gtest.h>


class BivariatePolynomialTests : public ::testing::Test
{
protected:
    BivariatePolynomial * Bp1, * Bp2;
    virtual void SetUp() {
        char p1[20];
        strcpy(p1, "4+3*x^2+2*x*y");
        char p2[20];
        strcpy(p2, "3+4*x^2*y^2+2*x*y");
        Bp1 = new BivariatePolynomial(p1, 2);
        Bp1->Print();
        Bp2 = new BivariatePolynomial(p2, 4);
        Bp2->Print();
    }
};



TEST_F(BivariatePolynomialTests, TestBackSubstitute1)
{
    double x = Bp1->backSubstitute(2, 1);
    EXPECT_EQ(20, x);
    x = Bp2->backSubstitute(3, 2);
    EXPECT_EQ(159, x);
}
