#include "ProblemSolver.h"
#include <lapacke.h>
#define NUM_OF_TRIES 3
using namespace std;
using namespace Eigen;

ProblemSolver::ProblemSolver(SylvesterPolynomial * SP, unsigned int B, unsigned int d):degree(d), mdIsInvertible(true), Solved(false)
{
    this->sylvesterPolynomial = SP;
    computeUpperBound(B);
    initMd(SP->getMd(this->degree), SP->getMatrixDimensions());
    //double k[] = {1,2,3,4,5,6,7,8,9,10};
    //Md = Map<MatrixXd>(k, 4, 4);
    //double * temp = *(SP->getMd(this->degree));
    //Md = Map<MatrixXd>(*(SP->getMd(this->degree)), SP->getMatrixDimensions(), SP->getMatrixDimensions());
    this->dimensionM = this->sylvesterPolynomial->getMatrixDimensions();
    cout<<"\n\n!!!! Md is !!!!\n"<<endl;
    cout << Md << endl<<endl;
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

}

bool ProblemSolver::isUnsolvable()
{

    srand(time(NULL));
    MatrixXd mI(this->dimensionM, this->dimensionM);
    for (int i = 0; i < NUM_OF_TRIES; ++i) {
        int randomY0 = rand() % 100 - 50;
        MatrixXd m(this->dimensionM, this->dimensionM);
        for (int j = 0; j < this->degree; ++j) {
            int powerY0 = randomY0;
            powerOf(powerY0, j);
            initMI(mI, this->sylvesterPolynomial->getMd(j), this->dimensionM);
            m += mI * powerY0;
        }
        powerOf(randomY0, this->degree);
        m += this->Md * randomY0;
        if(m.determinant() == 0)
            return true;
//        cout<<"m is : "<<endl<<m<<endl<<endl;
//        cout<<"determinant is : "<<m.determinant()<<endl;
    }
    return false;

}

void ProblemSolver::powerOf(int &  Num, int power)
{
    if(power == 0)
        Num = 1;
    for (int i = 0; i < power; ++i) {
        Num *= Num;
    }
}

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

bool ProblemSolver::solveStandardEigenProblem()
{

    if(degree < 1)
        return false;
    computeCompanion();
    cout<<this->Companion<<endl<<endl;

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
    for (int i = 0; i < eivals.rows(); ++i) {
        this->Solutions[i] = new Solution(eivecs(eivecs.rows() - 2,i).real()/eivecs(eivecs.rows() - 1,i).real(),
                                          -eivals(i).real(), multiplicity[i]);
        this->Solutions[i]->PrintSolution();
    }
    return true;
}
bool ProblemSolver::StandardNoSolution(MatrixXcd Eivecs, MatrixXcd Eivals, int i)
{
    if(Eivecs(Eivecs.rows() - 1, i).real() == 0 && Eivals(Eivecs.rows() - 1, i).imag() == 0)                  //no solution
        return true;
    if(isCloseToZero(Eivals(i).imag()))
        return false;
}

bool ProblemSolver::isCloseToZero(double Num)
{
    if(Num < 0.001)
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
        this->Companion.block(dimensionM*i,dimensionM*(i+1), dimensionM, dimensionM) = identity;
        for (int j = i + 2; j < degree ; ++j) {
            this->Companion.block(dimensionM*i,dimensionM*j, dimensionM, dimensionM) = zeroes;
        }
    }
    for (int i = 0; i < degree; ++i) {
        MatrixXd mI(dimensionM, dimensionM);
        this->initMI(mI, this->sylvesterPolynomial->getMd(i), dimensionM);
        this->Companion.block(dimensionM*(degree - 1), dimensionM*i, dimensionM, dimensionM) = - (this->Md.inverse() * mI);
    }
    return true;
}


bool ProblemSolver::solveGeneralizedEigenProblem()
{
    MatrixXd L0(this->degree * dimensionM, this->degree * dimensionM);
    MatrixXd L1(this->degree * dimensionM, this->degree * dimensionM);
    this->computeL0(L0);
    this->computeL1(L1);
/*    cout<<"----L0---"<<degree<<dimensionM<<endl;
    cout<<L0<<endl<<endl;
    cout<<"----L1----"<<endl;
    cout<<L1<<endl<<endl;*/
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
    cout<<"EigenValues are : "<<endl<<Eivals<<endl<<endl;
    cout<<"EigenVectors are : "<<endl<<Eivecs<<endl<<endl;
    cout<<"Roots"<<endl;
    cout<<"--------"<<endl;
    this->Solutions = new Solution * [Eivals.rows()];
    for (int i = 0; i < Eivals.rows(); ++i) {
        this->Solutions[i] = new Solution(Eivecs(this->degree * dimensionM - 2,i)/Eivecs(degree * dimensionM - 1,i),
                                          -Eivals(i, 0), multiplicity[i]);
        this->Solutions[i]->PrintSolution();
    }

    return false;
}

bool ProblemSolver::GeneralizedNoSolution(MatrixXd Eivecs, MatrixXd Eivals, int i)
{
    if(Eivecs(Eivecs.rows() - 1, i) == 0 && Eivecs(Eivecs.rows() - 1, i) == 0)                  //no solution
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

unsigned int ProblemSolver::computeUpperBound(unsigned int power)
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
    for (int j = 0; j < Eivals.rows(); ++j) {
        for (int l = 0; l < Eivals.rows(); ++l) {
            if((Eivals(j) == Eivals(l) && l != j) || this->StandardNoSolution(Eivecs, Eivals, l)) {
                if(!this->StandardNoSolution(Eivecs, Eivals, l))
                {
//                    cout<<"MULTIPLICITY ++ "<<endl;
//                    cout<<Eivecs(l)<<Eivals(l);
                    Multiplicity[j]++;
                }
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
    return ret;
}

bool ProblemSolver::removeSolsWithMultiplicityGeneralized(Eigen::MatrixXd & Eivecs, Eigen::MatrixXd & Eivals, int * Multiplicity)
{
    bool ret = false;
    int numRows = Eivals.rows() - 1;
    int numCols = Eivecs.cols() - 1;
    for (int j = 0; j < Eivals.rows(); ++j) {
        for (int l = 0; l < Eivals.rows(); ++l) {
            if((Eivals(j) == Eivals(l) && l != j ) || this->GeneralizedNoSolution(Eivecs, Eivals, l)) {
                if(!this->GeneralizedNoSolution(Eivecs, Eivals, l))
                {
//                    cout<<"MULTIPLICITY ++ "<<endl;
//                    cout<<Eivecs(l)<<Eivals(l,0);
                    Multiplicity[j]++;
                }
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
    return ret;
}