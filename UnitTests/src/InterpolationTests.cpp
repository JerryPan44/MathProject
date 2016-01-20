#include "../../Interpolation.h"
#include <gtest/gtest.h>
using namespace std;
class InterpolationTests : public ::testing::Test
{
protected:
    Interpolation* Inp;
    virtual void SetUp() {
        MatrixXd B(5,2);
	B<<-1, 0, 4, -1, -1, -5, -4, 2, 4, -4;
	cout << B << endl;
        Inp= new Interpolation(2,B);
    }
};

TEST_F(InterpolationTests, Testentry1)
{
    double x = Inp->entry(1,1,2);
    EXPECT_EQ(4, x);
    x = Inp->entry(4,2,2);
    EXPECT_EQ(256, x);
}

TEST_F(InterpolationTests, TestpowerOf1)
{
    double x = Inp->powerOf(6.2,2);
    EXPECT_EQ(38.44, x);
    x = Inp->powerOf(11.6,1);
    EXPECT_EQ(11.6, x);
}

TEST_F(InterpolationTests, TestcomputeM1)
{
    MatrixXd M = Inp->ComputeM();
    EXPECT_EQ(M(0,0), 1);
    EXPECT_EQ(M(1,3), 16);
    EXPECT_EQ(M(0,2), 0);
    EXPECT_EQ(M(0,5), 0);
    EXPECT_EQ(M(2,5), 25);
    EXPECT_EQ(M(2,1), -1);
    EXPECT_EQ(M(3,1), -4);
    EXPECT_EQ(M(4,4), -16);
    EXPECT_EQ(M(3,0), 1);
    EXPECT_EQ(M(2,0), 1);
}

TEST_F(InterpolationTests, Testfind1)
{
    ofstream equationsTxt;
    equationsTxt.open("InterpolationEquations.txt");
    BivariatePolynomial* a=Inp->find(equationsTxt);
    EXPECT_EQ(Inp->getRank(), 5);
    EXPECT_EQ(Inp->getK(), 5);
}

