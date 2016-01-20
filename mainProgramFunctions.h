#include <cstdio>
#include <cstring>
#include <unistd.h>
#include "eigen/Eigen/Dense"
#include "eigen/Eigen/Core"
#include "eigen/Eigen/SVD"
#include "Parse.h"
#include "BivariatePolynomial.h"
#include "SylvesterMatrix.h"
#include "Polynomial.h"
#include "SylvesterPolynomial.h"
#include "Matrix.h"
#include "ProblemSolver.h"
#include "ChangeOfVariableCoefficients.h"
#include "Interpolation.h"
#include "ErrorMargin.h"
#include <fstream>
#include <cmath>
using namespace std;
unsigned int getHiddenMaxDeg(SylvesterMatrix * SM, BivariatePolynomial * Bp1, BivariatePolynomial * Bp2);
void findFullSolutions(BivariatePolynomial * Bp1, BivariatePolynomial * Bp2, ProblemSolver * PS);
void cleanResources(ProblemSolver * PS, SylvesterMatrix * , BivariatePolynomial * , BivariatePolynomial * , SylvesterPolynomial * , MyMatrix *   );
void solveProblem(char * filename, int d1, int d2,int B, bool fromPoints);
void solveGeneratedProblem(int d1, int d2, int B);
void changeOfVariable(BivariatePolynomial * Bp1, BivariatePolynomial * Bp2,SylvesterMatrix * SM, ProblemSolver * PS, int B);
bool backSubstituteSols(BivariatePolynomial * Bp1, BivariatePolynomial * Bp2, ProblemSolver * PS, ofstream &);
bool backSubstituteSols(BivariatePolynomial * Bp1, BivariatePolynomial * Bp2, ProblemSolver * PS, ProblemSolver * PSNew);
