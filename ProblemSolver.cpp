#include "ProblemSolver.h"
#include <lapacke.h>
#include <cmath>
#define NUM_OF_TRIES 3
using namespace std;
using namespace Eigen;

ProblemSolver::ProblemSolver(SylvesterPolynomial * SP, unsigned int B, unsigned int d):degree(d), mdIsInvertible(true), Solved(false)
{
    this->sylvesterPolynomial = SP;
    computeUpperBound(B);
    initMd(SP->getMd(this->degree), SP->getMatrixDimensions());//Compute Md
    //double k[] = {1,2,3,4,5,6,7,8,9,10};
    //Md = Map<MatrixXd>(k, 4, 4);
    //double * temp = *(SP->getMd(this->degree));
    //Md = Map<MatrixXd>(*(SP->getMd(this->degree)), SP->getMatrixDimensions(), SP->getMatrixDimensions());
    this->dimensionM = this->sylvesterPolynomial->getMatrixDimensions();
    cout<<"\n\n!!!! Md is !!!!\n"<<endl;
    cout << Md << endl<<endl;			//1)
    this->k = this->computeStateIndicator();
    if(this->mdIsInvertible == true &&  k < upperBound )
        cout<<"k = "<<k<<" < Bound: non-singular Μd, standard eigenproblem"<<endl<<endl;
    if(this->mdIsInvertible == true && k > upperBound)
        cout<<"k = "<<k<<" > Bound: ill-conditioned Μd, generalized eigenproblem"<<endl<<endl;
    if(this->mdIsInvertible == false)
        cout<<"k --> infinity"<<" > Bound: ill-conditioned Μd, generalized eigenproblem"<<endl;

}

ProblemSolver::~ProblemSolver()
{
    for (int i = 0; i < this->numOfSolutions; ++i) {
        delete this->Solutions[i];
    }
}
// get determinant of the sylvester Matrix for a random Y0 3 times if 0 for all the Y0 means that probably the system has no solution
bool ProblemSolver::isUnsolvable()
{									//5)[BONUS]

    srand(time(NULL));
    MatrixXd mI(this->dimensionM, this->dimensionM);
    for (int i = 0; i < NUM_OF_TRIES; ++i) {
        int randomY0 = rand() % 100 - 50;
        MatrixXd m = MatrixXd::Zero(this->dimensionM, this->dimensionM);
        for (int j = 0; j < this->degree; ++j) {
            int powerY0 = randomY0;
            powerOf(powerY0, j);
            initMI(mI, this->sylvesterPolynomial->getMd(j), this->dimensionM);
            m += mI * powerY0;
        }
        powerOf(randomY0, this->degree);
        m += this->Md * randomY0;
        if(m.determinant() != 0)
            return false;
//        cout<<"m is : "<<endl<<m<<endl<<endl;
//        cout<<"determinant is : "<<m.determinant()<<endl;
    }
    return true;
}

void ProblemSolver::powerOf(int &  Num, int power)
{
    if(power == 0)
        Num = 1;
    for (int i = 0; i < power; ++i) {
        Num *= Num;
    }
}
//solve the problem
bool ProblemSolver::Solve()
{
    if(this->isUnsolvable())
        return false;
    if(Solved == true)
        return false;
    Solved = true;
    if( k < upperBound && this->mdIsInvertible == true)
        return this->solveStandardEigenProblem();
    else
        return this->solveGeneralizedEigenProblem();
}

bool ProblemSolver::solveStandardEigenProblem()		//2)
{

    if(degree < 1)
        return false;
    computeCompanion();

    EigenSolver<MatrixXd> es(this->Companion);

    MatrixXcd eivecs = es.eigenvectors();
    MatrixXcd eivals = es.eigenvalues();
//    cout<<"EigenValues are : "<<endl<<eivals<<endl<<endl;
//    cout<<"EigenVectors are : "<<endl<<eivecs<<endl<<endl;
    int * multiplicity = new int[eivals.rows()];
    for (int j = 0; j < eivals.rows(); ++j) {
        multiplicity[j] = 1;
    }
    removeSolsWithMultiplicityStandard(eivecs, eivals, multiplicity);
    cout<<"Roots"<<endl;
    cout<<"----------"<<endl;
    this->Solutions = new Solution*[eivals.rows()];
    this->numOfSolutions = eivals.rows();
    for (int i = 0; i < eivals.rows(); ++i) {
        this->Solutions[i] = new Solution(eivecs(eivecs.rows() - 2,i).real()/eivecs(eivecs.rows() - 1,i).real(),
                                          eivals(i).real(), multiplicity[i]);
        this->Solutions[i]->PrintSolution();
    }
    return true;
}
//no solution if imaginary part is not zero and eivec cannot be normalized
bool ProblemSolver::StandardNoSolution(MatrixXcd Eivecs, MatrixXcd Eivals, int i, int * Multiplicity)
{
    if(Eivecs(Eivecs.rows() - 1, i).real() == 0 && Eivals(Eivecs.rows() - 1, i).imag() == 0 && Multiplicity[i] == 1)                  //no solution
        return true;
    if(isCloseToZero(Eivals(i).imag()))
        return false;
    return true;
}

bool ProblemSolver::isCloseToZero(double Num)
{
    if(abs(Num) < 0.001)
        return true;
    return false;
}


bool ProblemSolver::computeCompanion()
{
    this->Companion = MatrixXd(this->degree * dimensionM, this->degree * dimensionM);
    MatrixXd identity = MatrixXd::Identity(dimensionM, dimensionM);
    MatrixXd zeroes = MatrixXd::Zero(dimensionM, dimensionM);
    for (int i = 0; i < degree - 1; ++i) {
        for (int j = 0; j < i+1; ++j) {
            this->Companion.block(dimensionM*i,dimensionM*j, dimensionM, dimensionM) = zeroes;
        }
        this->Companion.block(dimensionM*i,dimensionM*(i+1), dimensionM, dimensionM) = identity;        //identity on the right of the diagonal 0 everywhere else
        for (int j = i + 2; j < degree ; ++j) {
            this->Companion.block(dimensionM*i,dimensionM*j, dimensionM, dimensionM) = zeroes;
        }
    }
    for (int i = 0; i < degree; ++i) {
        MatrixXd mI(dimensionM, dimensionM);
        this->initMI(mI, this->sylvesterPolynomial->getMd(i), dimensionM);          //get ith degree matrix
        this->Companion.block(dimensionM*(degree - 1), dimensionM*i, dimensionM, dimensionM) = - (this->Md.inverse() * mI);         //last row of companion matrix
    }
    return true;
}


bool ProblemSolver::solveGeneralizedEigenProblem()		//3)
{
    MatrixXd L0(this->degree * dimensionM, this->degree * dimensionM);
    MatrixXd L1(this->degree * dimensionM, this->degree * dimensionM);
    this->computeL0(L0);
    this->computeL1(L1);
//    cout<<"----L0---"<<degree<<dimensionM<<endl;
//    cout<<L0<<endl<<endl;
//    cout<<"----L1----"<<endl;
//    cout<<L1<<endl<<endl;
    MatrixXd Eivecs(this->degree * dimensionM, this->degree * dimensionM);
    MatrixXd Eivals(this->degree * dimensionM, 3);
    LapackeSolveGeneralizedEigenProblem(L0, L1, Eivecs, Eivals);

//    cout<<"EigenValues are : "<<endl<<Eivals<<endl<<endl;
//    cout<<"EigenVectors are : "<<endl<<Eivecs<<endl<<endl;
    int * multiplicity = new int[this->degree * dimensionM];
    for (int j = 0; j < this->degree * dimensionM; ++j) {
        multiplicity[j] = 1;
    }
    removeSolsWithMultiplicityGeneralized(Eivecs, Eivals, multiplicity);
//    cout<<"EigenValues are : "<<endl<<Eivals<<endl<<endl;
//    cout<<"EigenVectors are : "<<endl<<Eivecs<<endl<<endl;

    cout<<"Roots"<<endl;
    cout<<"--------"<<endl;
    this->Solutions = new Solution * [Eivals.rows()];
    this->numOfSolutions = Eivals.rows();
    for (int i = 0; i < Eivals.rows(); ++i) {
        this->Solutions[i] = new Solution(Eivecs(this->degree * dimensionM - 2,i)/Eivecs(degree * dimensionM - 1,i),
                                          -(Eivals(i, 0)/Eivals(i, 2)), multiplicity[i]);
        this->Solutions[i]->PrintSolution();
    }
    return false;
}
/// if last pos of eivalue has 0 denominator or eivec cannot be normalized 0 on last pos return true else return false
bool ProblemSolver::GeneralizedNoSolution(MatrixXd Eivecs, MatrixXd Eivals, int i, int * Multiplicity)
{
    if(Eivals(i, 2) == 0 || abs(Eivals(i,2)) < 0.00001)
        return true;
    if(Eivecs(Eivecs.rows() - 1, i) == 0 && Multiplicity[i] == 1)                  //no solution
        return true;
    if(isCloseToZero(Eivals(i, 1)))       //if the eigenvalue has no imaginary part
        return false;
    return true;
}

bool  ProblemSolver::LapackeSolveGeneralizedEigenProblem(MatrixXd& A, MatrixXd & B, MatrixXd& Eivecs, MatrixXd& Eivals)
{
    int N = A.cols();
    if(B.cols() != N ||   A.rows() != N || B.rows() != N)
        return false;
    int LDA = A.outerStride(); int LDB =  B.outerStride(); int LDV = Eivecs.outerStride(); int INFO = 0;

    double * alphar = Eivals.col(0).data();
    double * alphai = Eivals.col(1).data();
    double * beta = Eivals.col(2).data();
    INFO = LAPACKE_dggev(LAPACK_COL_MAJOR, 'N', 'V', N, A.data(), LDA, B.data(), LDB, alphar, alphai, beta, 0, LDV, Eivecs.data(), LDV);
    return INFO == 0;

}


bool ProblemSolver::computeL0(MatrixXd & L0)
{
    MatrixXd identity = - MatrixXd::Identity(dimensionM, dimensionM);
    MatrixXd zeroes = MatrixXd::Zero(dimensionM, dimensionM);
    for (int i = 0; i < degree - 1; ++i) {                          //for every row except the last
        L0.block(dimensionM*i,dimensionM*(i+1), dimensionM, dimensionM) = identity;         //identity is on the block matrix after the diagonal
        for (int j = 0; j < i + 1; ++j) {
            L0.block(dimensionM*i,dimensionM*j, dimensionM, dimensionM) = zeroes;       //set to  zero   before the identity
        }
        for (int j = i + 2; j < degree; ++j) {
            L0.block(dimensionM*i,dimensionM*j, dimensionM, dimensionM) = zeroes;   //set to zero  after the diagonal
        }
    }
    for (int i = 0; i < degree; ++i) {
        MatrixXd mI(dimensionM, dimensionM);
        initMI(mI, this->sylvesterPolynomial->getMd(i), dimensionM);
        L0.block(dimensionM*(degree - 1), dimensionM*i, dimensionM, dimensionM) = mI;
    }
    return true;
}

bool ProblemSolver::computeL1(MatrixXd & L1)
{
    const int dimensionM = this->sylvesterPolynomial->getMatrixDimensions();
    MatrixXd identity = MatrixXd::Identity(dimensionM, dimensionM);
    MatrixXd zeroes = MatrixXd::Zero(dimensionM, dimensionM);
    for (int i = 0; i < degree; ++i) {                          //for every row except the last
        if(i == degree - 1)
            L1.block(dimensionM*i,dimensionM*i, dimensionM, dimensionM) = this->Md;         //if on last row diagonal is Md
        else
            L1.block(dimensionM*i,dimensionM*i, dimensionM, dimensionM) = identity;         //identity is on the block matrix diagonal
        for (int j = 0; j < i; ++j) {
            L1.block(dimensionM*i,dimensionM*j, dimensionM, dimensionM) = zeroes;       //set to  zero   before the identity
        }
        for (int j = i + 1; j < degree ; ++j) {
            L1.block(dimensionM*i,dimensionM*j, dimensionM, dimensionM) = zeroes;   //set to zero  after the diagonal
        }
    }
    return true;
}
//initialize Md
bool ProblemSolver::initMd(double ** tempMD, int matrixDimensions)
{
    if(tempMD == NULL)
        return false;
    this->Md = MatrixXd(matrixDimensions, matrixDimensions);
    for (int i = 0; i < matrixDimensions; ++i) {
        for (int j = matrixDimensions - 1; j >= 0 ; --j) {
            this->Md(i, matrixDimensions - 1 - j) = tempMD[i][j];
        }
    }
    return true;
}

bool ProblemSolver::initMI(MatrixXd& mI, double ** tempMD, int matrixDimensions)
{
    if(tempMD == NULL)
        return false;
    for (int i = 0; i < matrixDimensions; ++i) {
        for (int j = matrixDimensions - 1; j >= 0 ; --j) {
            mI(i, matrixDimensions - 1 - j) = tempMD[i][j];
        }
    }
    return true;
}
//get max/min singular value if min is 0 then isInvertible = false
double ProblemSolver::computeStateIndicator()
{
    JacobiSVD<MatrixXd> svd(Md, ComputeFullU);
    JacobiSVD<MatrixXd>::SingularValuesType svds = svd.singularValues();
//    cout<<svds<<endl;
    double max = svds(0);
    double min = svds(0);
    for (int i = 1; i < svds.rows(); ++i) {
        if(svds(i) > max)
            max = svds(i);
        if(svds(i) < min)
            min = svds(i);
    }
    if(min == 0) {
            mdIsInvertible = false;
        return 1;
    }

    return max/min;
    // return eivals[0];

}

unsigned int ProblemSolver::computeUpperBound(unsigned int power)//Bound B
{
    this->upperBound = 10;
    for (int i = 0; i < power - 1; ++i) {
        this->upperBound *= 10;
    }
}
bool ProblemSolver::removeSolsWithMultiplicityStandard(Eigen::MatrixXcd & Eivecs, Eigen::MatrixXcd & Eivals, int * Multiplicity)
{
    bool ret = false;
    int numRows = Eivals.rows() - 1;
    int numCols = Eivecs.cols() - 1;
    for (int j = 0; j < Eivals.rows(); ++j) {           //for each eigen value if it is repeated resize the eivecs and eivals matrix and remove the duplicates
        for (int l = 0; l < Eivals.rows(); ++l) {
            if((Eivals(j).real() == Eivals(l).real() && Eivals(j).imag() == Eivals(l).imag() && l != j)) {
                Multiplicity[j]++;
                if( l < numRows )
                    Eivals.block(l, 0, numRows - l, Eivals.cols()) = Eivals.block(l + 1, 0, numRows - l, Eivals.cols());
                if( l < numCols )
                    Eivecs.block(0,l, Eivecs.rows(), numCols - l) = Eivecs.block(0,l+1,Eivecs.rows(),numCols-l);
                Eivals.conservativeResize(numRows,Eivals.cols());
                Eivecs.conservativeResize(Eivecs.rows(),numCols);
                numRows = Eivals.rows() - 1;
                numCols = Eivecs.cols() - 1;

                ret = true;
            }
        }    
    }
    for (int j = 0; j < Eivals.rows(); ++j) {
        if(this->StandardNoSolution(Eivecs, Eivals, j, Multiplicity))           //if the eivector and eigenvalues do not gratify the constraints we have set for them remove the solution
        {
            if( j < numRows )
                Eivals.block(j, 0, numRows - j, Eivals.cols()) = Eivals.block(j + 1, 0, numRows - j, Eivals.cols());
            if( j < numCols )
                Eivecs.block(0,j, Eivecs.rows(), numCols - j) = Eivecs.block(0,j+1,Eivecs.rows(),numCols-j);
            Eivals.conservativeResize(numRows,Eivals.cols());
            Eivecs.conservativeResize(Eivecs.rows(),numCols);
            numRows = Eivals.rows() - 1;
            numCols = Eivecs.cols() - 1;
            continue;
        }
    }
    return ret;
}
///Same constraint logic as the standard eigen problem
bool ProblemSolver::removeSolsWithMultiplicityGeneralized(Eigen::MatrixXd & Eivecs, Eigen::MatrixXd & Eivals, int * Multiplicity)
{
    bool ret = false;
    int numRows = Eivals.rows() - 1;
    int numCols = Eivecs.cols() - 1;
    for (int j = 0; j < Eivals.rows(); ++j) {
        for (int l = 0; l < Eivals.rows(); ++l) {
            if(Eivals(j,0) == Eivals(l,0) && Eivals(j,2) == Eivals(l,2) && l != j ) {
                Multiplicity[j]++;
//                if(this->GeneralizedNoSolution(Eivecs, Eivals, l))
//                    cout<<Eivals(j, 0)<<"/"<<Eivals(j,2)<<"  "<<Eivals(l,0)<<"/"<<Eivals(l,2)<<endl;
                if( l < numRows )
                    Eivals.block(l, 0, numRows - l, Eivals.cols()) = Eivals.block(l + 1, 0, numRows - l, Eivals.cols());
                if( l < numCols )
                    Eivecs.block(0,l, Eivecs.rows(), numCols - l) = Eivecs.block(0,l+1,Eivecs.rows(),numCols-l);
                Eivals.conservativeResize(numRows,Eivals.cols());
                Eivecs.conservativeResize(Eivecs.rows(),numCols);
                numRows = Eivals.rows() - 1;
                numCols = Eivecs.cols() - 1;
                ret = true;
            }
        }
    }
    for (int j = 0; j < Eivals.rows(); ++j) {

        if(this->GeneralizedNoSolution(Eivecs, Eivals, j, Multiplicity))
        {
            if( j < numRows )
                Eivals.block(j, 0, numRows - j, Eivals.cols()) = Eivals.block(j + 1, 0, numRows - j, Eivals.cols());
            if( j < numCols )
                Eivecs.block(0,j, Eivecs.rows(), numCols - j) = Eivecs.block(0,j+1,Eivecs.rows(),numCols-j);
            Eivals.conservativeResize(numRows,Eivals.cols());
            Eivecs.conservativeResize(Eivecs.rows(),numCols);
            numRows = Eivals.rows() - 1;
            numCols = Eivecs.cols() - 1;
            continue;
        }
    }
    return ret;
}
///get the old y back by replacing with z (y = t1z+t2/t3x+t4)
bool ProblemSolver::substituteChangeOfVariable(ChangeOfVariableCoefficients * coefs)			//4)
{
    for (int j = 0; j < this->numOfSolutions; ++j) {
        double y = this->Solutions[j]->getY();
        double x = this->Solutions[j]->getX();
        int mul = this->Solutions[j]->getMultiplicity();
//        y = (coefs->t2 - y * coefs->t4)/(y * coefs->t3 - coefs->t1);
        y = (coefs->t1 * y + coefs->t2) / (coefs->t3 * y + coefs->t4);
        delete this->Solutions[j];
        this->Solutions[j] = new Solution(x, y, mul);
        this->Solutions[j]->PrintSolution();
    }
}
