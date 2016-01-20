
#include <iostream>
#include <gtest/gtest.h>
#include "InterpolationTests.cpp"
#include "PolynomialTests.cpp"
#include "SylvesterMatrixTests.cpp"
#include "ProblemSolverTests.cpp"
#include "BivariatePolynomialTests.cpp"
#include "MainProgramTests.cpp"
int main(int argc, char * argv[])
{

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}