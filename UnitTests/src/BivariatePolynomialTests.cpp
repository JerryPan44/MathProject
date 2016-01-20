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
//        Bp1->Print();
        Bp2 = new BivariatePolynomial(p2, 4);
//        Bp2->Print();
    }
};



TEST_F(BivariatePolynomialTests, TestBackSubstitute1)
{
    double x = Bp1->backSubstituteXandY(2, 1);
    EXPECT_EQ(20, x);
    x = Bp2->backSubstituteXandY(3, 2);
    EXPECT_EQ(159, x);
}

TEST_F(BivariatePolynomialTests, TestBackSubstitute2)
{
    char p1[20];
    strcpy(p1, "1+y*x-x+y");
    char p2[20];
    strcpy(p2, "1+y*x^2-x^2+2*y*x+y+x");
    Bp1 = new BivariatePolynomial(p1, 2);
    Bp2 = new BivariatePolynomial(p2, 3);
    Bp1->Print();
    double x = Bp1->backSubstituteXandY(0, -1);
    EXPECT_EQ(0, x);
    x = Bp2->backSubstituteXandY(0, -1);
    EXPECT_EQ(0, x);
}
TEST_F(BivariatePolynomialTests, TestBackSubstitute3)
{
  	int i=2;
	char a[]="2.73333333+2.25*x+0.50*y-4.8333333*x^2+0.100*y^2";
	BivariatePolynomial* Bp1 = new BivariatePolynomial(a,i);
	Bp1->Print();
}