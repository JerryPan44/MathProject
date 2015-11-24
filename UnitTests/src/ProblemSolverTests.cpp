#include "../../SylvesterMatrix.h"
#include "../../Polynomial.h"
#include "../../BivariatePolynomial.h"
#include "../../ProblemSolver.h"
#include "../../Matrix.h"
#include <gtest/gtest.h>

using namespace Eigen;
using namespace std;
class ProblemSolverTests : public ::testing::Test {
protected:
    SylvesterMatrix * SM;
    SylvesterPolynomial * SP, * SP2;
    BivariatePolynomial * Bp1, * Bp2;
    ProblemSolver * PS, * PS2;
    virtual void SetUp() {


    }


};

TEST_F(ProblemSolverTests, testMd)
{
    char p1[20];
    strcpy(p1, "4+3*x^2+2*x*y");
    char p2[20];
    strcpy(p2, "3+4*x^2+5*x*y");
    Bp1 = new BivariatePolynomial(p1, 2);
    //        Bp1->Print();
    Bp2 = new BivariatePolynomial(p2, 2);
    //        Bp2->Print();
    SM = new SylvesterMatrix(Bp1, Bp2);
    SP = new SylvesterPolynomial(SM->getHiddenDeg(), SM->getColDimension());
    SP->SMatrixToSPolynomial(SM);
    //        SP->Print();
    PS = new ProblemSolver(SP, 10, SP->getDegree());
    double ** m = SP->getMd(SP->getDegree());
    MatrixXd m2 = PS->getMd();
    for (int i = 0; i < SP->getMatrixDimensions(); ++i) {
        for (int j = 0; j < SP->getMatrixDimensions(); ++j) {
            EXPECT_EQ(m[i][SP->getMatrixDimensions() - 1 - j], m2(i, j));
        }
    }
}

TEST_F(ProblemSolverTests, testMd2)
{
    char p2[20];
    char p1[20];
    strcpy(p1, "-3+2*y*x+y");
    strcpy(p2, "5+y*x^2+4*x-y");
    Bp1 = new BivariatePolynomial(p1, 2);
    Bp2 = new BivariatePolynomial(p2, 2);
    SM = new SylvesterMatrix(Bp1, Bp2);
    SP2 = new SylvesterPolynomial(SM->getHiddenDeg(), SM->getColDimension());
    SP2->SMatrixToSPolynomial(SM);
    PS2 = new ProblemSolver(SP2, 10, SP2->getDegree());
    MatrixXd m = PS2->getMd();
    EXPECT_EQ(2, m(0,0));
    EXPECT_EQ(1, m(0,1));
    EXPECT_EQ(2, m(1,1));
    EXPECT_EQ(1, m(1,2));
    EXPECT_EQ(1, m(2,0));
    EXPECT_EQ(-1, m(2,2));
    int k = PS2->getStateIndicator();
    EXPECT_EQ(5, k);
    PS2->Solve();
}

TEST_F(ProblemSolverTests, testMd3)
{
    char p1[20], p2[20];
    strcpy(p1, "1+y*x-x+y");
    strcpy(p2, "1+y*x^2-x^2+2*y*x+y+x");
    Bp1 = new BivariatePolynomial(p1, 2);
    Bp2 = new BivariatePolynomial(p2, 3);
    SM = new SylvesterMatrix(Bp1, Bp2);
    SP2 = new SylvesterPolynomial(SM->getHiddenDeg(), SM->getRowDimension());
    SP2->SMatrixToSPolynomial(SM);
    PS2 = new ProblemSolver(SP2, 7, SP2->getDegree());
    MatrixXd m = PS2->getMd();
    EXPECT_EQ(1, m(0,0));
    EXPECT_EQ(1, m(0,1));
    EXPECT_EQ(1, m(1,1));
    EXPECT_EQ(1, m(1,2));
    EXPECT_EQ(1, m(2,0));
    EXPECT_EQ(2, m(2,1));
    EXPECT_EQ(1, m(2,2));
    unsigned int k = PS2->getStateIndicator();
    EXPECT_EQ(4265241472, k);
//    SP2->Print();
    PS2->Solve();
}

TEST_F(ProblemSolverTests, testSolution)
{
    char p1[20], p2[20];
    strcpy(p1, "1+x^2-y");
    strcpy(p2, "-8+x^2+5*y");
    Bp1 = new BivariatePolynomial(p1, 2);
    Bp2 = new BivariatePolynomial(p2, 2);
    SM = new SylvesterMatrix(Bp1, Bp2);
    SP2 = new SylvesterPolynomial(SM->getHiddenDeg(), SM->getRowDimension());
    SP2->SMatrixToSPolynomial(SM);
    PS2 = new ProblemSolver(SP2, 7, SP2->getDegree());

    MatrixXd m = PS2->getMd();
    EXPECT_EQ(-1, m(0,2));
    EXPECT_EQ(-1, m(1,3));
    EXPECT_EQ(5, m(2,2));
    EXPECT_EQ(5, m(3,3));
    unsigned int k = PS2->getStateIndicator();
//    EXPECT_EQ(4265241472, k);
//    SP2->Print();
    PS2->Solve();
}

