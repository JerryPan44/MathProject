#ifndef _PROBLEM_SOLVER_
#define _PROBLEM_SOLVER_

#include <iostream>
#include "SylvesterPolynomial.h"
#include "eigen/Eigen/Dense"
#include "eigen/Eigen/Core"
#include "eigen/Eigen/SVD"
#include "eigen/Eigen/Eigenvalues"
#include "Solution.h"
#include "ChangeOfVariableCoefficients.h"
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MyMatrixXd;
class ProblemSolver
{
    bool Solved;
    Solution ** Solutions;
    int numOfSolutions;
    Eigen::MatrixXd Companion;
    Eigen::MatrixXd Md;
    unsigned int upperBound;
    double k;
    int degree, dimensionM;
    bool mdIsInvertible;
    bool solveStandardEigenProblem();
    bool solveGeneralizedEigenProblem();
    SylvesterPolynomial * sylvesterPolynomial;
    bool computeL0(Eigen::MatrixXd & L0);
    bool computeL1(Eigen::MatrixXd & L1);
    bool computeCompanion();
    bool initMd(double ** tempMD, int matrixDimensions);
    unsigned int computeUpperBound(unsigned int power);
    double computeStateIndicator();
    bool initMI(Eigen::MatrixXd& mI, double ** tempMD, int matrixDimensions);
    bool LapackeSolveGeneralizedEigenProblem(Eigen::MatrixXd& A, Eigen::MatrixXd & B, Eigen::MatrixXd& Eivecs, Eigen::MatrixXd& Eivals);
    bool isUnsolvable();
    bool GeneralizedNoSolution(Eigen::MatrixXd Eivecs, Eigen::MatrixXd Eivals, int i, int *);
    bool StandardNoSolution(Eigen::MatrixXcd & Eivecs, Eigen::MatrixXcd & Eivals, int i, int *);
    void powerOf(int& Num, int power);
    bool isCloseToZero(double Num);
    bool removeSolsWithMultiplicityStandard(Eigen::MatrixXcd &, Eigen::MatrixXcd &, int * );
    bool removeSolsWithMultiplicityGeneralized(Eigen::MatrixXd &, Eigen::MatrixXd &, int * );
public:

    bool getMdIsInvertible()
    {
        return this->mdIsInvertible;
    }

    double getStateIndicator()
    {
        return k;
    }
    Eigen::MatrixXd getMd()
    {
        return Md;
    }

    ProblemSolver(SylvesterPolynomial *, unsigned int, unsigned int);
    ~ProblemSolver();
    bool Solve();
    int getNumOfSols()
    {
        return this->numOfSolutions;
    }

    bool substituteChangeOfVariable(ChangeOfVariableCoefficients *);
    Solution * getSolution(int i)
    {
        if(this->numOfSolutions - 1 < i || i < 0)
            return  NULL;
        return this->Solutions[i];
    }



};


#endif