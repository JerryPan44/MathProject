#include "ProblemSolver.h"
#include <lapacke.h>
#define NUM_OF_TRIES 3
using namespace std;
using namespace Eigen;

ProblemSolver::ProblemSolver(SylvesterPolynomial * SP, unsigned int B, unsigned int d):degree(d), mdIsInvertible(true), Solved(false)
{														//Find Md type (non-singular,ill-consitioned) and if gen. or std eigenproblem
    this->sylvesterPolynomial = SP;
    computeUpperBound(B);											//compute 10^B (Upper Bound)
    initMd(SP->getMd(this->degree), SP->getMatrixDimensions());							//Compute Md
    this->dimensionM = this->sylvesterPolynomial->getMatrixDimensions();
    cout<<"\n\n!!!! Md is !!!!\n"<<endl;
    cout << Md << endl<<endl;											//Print Md
    this->k = this->computeStateIndicator();									//Compute k
    if(this->mdIsInvertible == true &&  k < upperBound )							//Is eigenproblem standard or generalised
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
{									//Judge if system can be solved and answer to 5)[BONUS] 
    
    srand(time(NULL));
    MatrixXd mI(this->dimensionM, this->dimensionM);
    for (int i = 0; i < NUM_OF_TRIES; ++i) {				//Reapeat up to 3 times (NUM_OF_TRIES=3)
        int randomY0 = rand() % 100 - 50;				//Pick a random number over {-50,...,50}
        MatrixXd m(this->dimensionM, this->dimensionM);
        for (int j = 0; j < this->degree; ++j) {
            int powerY0 = randomY0;
            powerOf(powerY0, j);					//(Y0)^j
            initMI(mI, this->sylvesterPolynomial->getMd(j), this->dimensionM);
            m += mI * powerY0;						//m= mI*y0^1 + mI*y0^2 +...+ mI*y0^(degree-1)
        }
        powerOf(randomY0, this->degree);				//(Y0)^(degree)
        m += this->Md * randomY0;
        if(m.determinant() == 0)					//det [m(Y0)]
            return true;						//Edw nomizw prepei na allaksei (des Allages-Protaseis)
//        cout<<"m is : "<<endl<<m<<endl<<endl;
//        cout<<"determinant is : "<<m.determinant()<<endl;
    }
    return false;

}

void ProblemSolver::powerOf(int &  Num, int power)			//Compute Num^(power)
{
    if(power == 0)
        Num = 1;
    for (int i = 0; i < power; ++i) {
        Num *= Num;
    }
}

bool ProblemSolver::Solve()
{									//Judge if system is: solved, unsolvable, a standard eigenproblem, a generalised eigenproblem
    if(this->isUnsolvable())						//System cannot be solved
        return false;
    if(Solved == true)							//System is solved
        return false;
    Solved = true;
    if( k < upperBound && this->mdIsInvertible == true)
        return this->solveStandardEigenProblem();			//Solve with the standard eigenproblem method
    else
        return this->solveGeneralizedEigenProblem();			//Solve with the generalised eigenproblem method
}

bool ProblemSolver::solveStandardEigenProblem()				//Solve with the standard eigenproblem method
{
    if(degree < 1)
        return false;
    computeCompanion();							//compute Companion matrix
    EigenSolver<MatrixXd> es(this->Companion);				//Give created Companion matrix to eigen lib
    MatrixXcd eivecs = es.eigenvectors();				//Eigenvectors
    MatrixXcd eivals = es.eigenvalues();				//Eigenvalues
//    cout<<"EigenValues are : "<<endl<<eivals<<endl<<endl;
//    cout<<"EigenVectors are : "<<endl<<eivecs<<endl<<endl;
    int * multiplicity = new int[eivals.rows()];
    for (int j = 0; j < eivals.rows(); ++j) {
        multiplicity[j] = 1;						//For every returned eigenvalue its multiplicity becomes 1
    }
    removeSolsWithMultiplicityStandard(eivecs, eivals, multiplicity);	//
    cout<<"Roots"<<endl;
    cout<<"----------"<<endl;
    this->Solutions = new Solution*[eivals.rows()];			//commit pointers in number equal to the number of solutions
    this->numOfSolutions = eivals.rows();				//Number of solutions updated to the class
    for (int i = 0; i < eivals.rows(); ++i) {
        this->Solutions[i] = new Solution(eivecs(eivecs.rows() - 2,i).real()/eivecs(eivecs.rows() - 1,i).real(),
                                          eivals(i).real(), multiplicity[i]);	//add solutions
        this->Solutions[i]->PrintSolution();					//and print them
    }
    return true;
}

bool ProblemSolver::StandardNoSolution(MatrixXcd Eivecs, MatrixXcd Eivals, int i, int * Multiplicity)
{												//judge if Generalised eigenproblem has no solution given i postition

    if(Eivecs(Eivecs.rows() - 1, i).real() == 0 && Eivals(Eivecs.rows() - 1, i).imag() == 0 && Multiplicity[i] == 1)
        return true;										//eigenvalue cannot be canonized
    if(isCloseToZero(Eivals(i).imag()))
        return false;										//if the eigenvalue has imaginary part < 10^-3
    return true;
}

bool ProblemSolver::isCloseToZero(double Num)								//A number is considered closed enough to zero if it is < 10^-3
{
    if(Num < 0.001)
        return true;
    return false;
}


bool ProblemSolver::computeCompanion()									//Create and compute the Companion Matrix
{
    this->Companion = MatrixXd(this->degree * dimensionM, this->degree * dimensionM);			//create Companion matrix
    MatrixXd identity = MatrixXd::Identity(dimensionM, dimensionM);					//indetity is the intetity matrix of dimM x dimM
    MatrixXd zeroes = MatrixXd::Zero(dimensionM, dimensionM);						//zeroes is the zero matrix of dimM x dimM
    for (int i = 0; i < degree - 1; ++i) {								//For every row except the last
        for (int j = 0; j < i+1; ++j) {
            this->Companion.block(dimensionM*i,dimensionM*j, dimensionM, dimensionM) = zeroes;		//first i blocks are zero matrices
        }
        this->Companion.block(dimensionM*i,dimensionM*(i+1), dimensionM, dimensionM) = identity;	//the (i+1)-th block is an identity matrix
        for (int j = i + 2; j < degree ; ++j) {
            this->Companion.block(dimensionM*i,dimensionM*j, dimensionM, dimensionM) = zeroes;		//the rest are zero matrices
        }
    }
    for (int i = 0; i < degree; ++i) {
        MatrixXd mI(dimensionM, dimensionM);
        this->initMI(mI, this->sylvesterPolynomial->getMd(i), dimensionM);				//Create Mi, i=0,1,...
        this->Companion.block(dimensionM*(degree - 1), dimensionM*i, dimensionM, dimensionM) = - (this->Md.inverse() * mI);
    }													//Last row of blocks has -Md^(-1)*Mi
    return true;
}


bool ProblemSolver::solveGeneralizedEigenProblem()					//Solve with the genralised eigenproblem method !!!(needs rev)!!!!
{
    MatrixXd L0(this->degree * dimensionM, this->degree * dimensionM);			//create L0 matrix
    MatrixXd L1(this->degree * dimensionM, this->degree * dimensionM);			//create L1 matrix
    this->computeL0(L0);								//compute L0 matrix
    this->computeL1(L1);								//compute L1 matrix
//    cout<<"----L0---"<<degree<<dimensionM<<endl;
//    cout<<L0<<endl<<endl;
//    cout<<"----L1----"<<endl;
//    cout<<L1<<endl<<endl;
    MatrixXd Eivecs(this->degree * dimensionM, this->degree * dimensionM);		//Eigenvectors
    MatrixXd Eivals(this->degree * dimensionM, 3);					//Eigenvalues
    LapackeSolveGeneralizedEigenProblem(L0, L1, Eivecs, Eivals);			//Solve the generalized eigenproblem using Lapacke library

//    cout<<"EigenValues are : "<<endl<<Eivals<<endl<<endl;
//    cout<<"EigenVectors are : "<<endl<<Eivecs<<endl<<endl;
    int * multiplicity = new int[this->degree * dimensionM];
    for (int j = 0; j < this->degree * dimensionM; ++j) {
        multiplicity[j] = 1;								//For every returned eigenvalue its multiplicity becomes 1
    }
    removeSolsWithMultiplicityGeneralized(Eivecs, Eivals, multiplicity);		//Remove solutions with multiplicity >1
//    cout<<"EigenValues are : "<<endl<<Eivals<<endl<<endl;
//    cout<<"EigenVectors are : "<<endl<<Eivecs<<endl<<endl;

    cout<<"Roots"<<endl;
    cout<<"--------"<<endl;
    this->Solutions = new Solution * [Eivals.rows()];					//commit pointers in number equal to the number of solutions
    this->numOfSolutions = Eivals.rows();						//Number of solutions updated to the class
    for (int i = 0; i < Eivals.rows(); ++i) {
        this->Solutions[i] = new Solution(Eivecs(this->degree * dimensionM - 2,i)/Eivecs(degree * dimensionM - 1,i),
                                          -(Eivals(i, 0)/Eivals(i, 2)), multiplicity[i]);	//add solutions
        this->Solutions[i]->PrintSolution();							//and print them
    }
    return false;
}

bool ProblemSolver::GeneralizedNoSolution(MatrixXd Eivecs, MatrixXd Eivals, int i, int * Multiplicity)
{											//judge if Generalised eigenproblem has no solution given i postition
    if(Eivals(i, 2) == 0 || Eivals(i,2) < 0.00001)					//eigenvalue is 0 or close to 0 (<10^-5)
        return true;
    if(Eivecs(Eivecs.rows() - 1, i) == 0 && Multiplicity[i] == 1)			//eigenvector cannot be canonized
        return true;
    if(isCloseToZero(Eivals(i, 1)))							//if the eigenvalue has imaginary part < 10^-3
        return false;
    return true;
}

bool  ProblemSolver::LapackeSolveGeneralizedEigenProblem(MatrixXd& A, MatrixXd & B, MatrixXd& Eivecs, MatrixXd& Eivals)
{											//Solve the generalized eigenproblem using Lapacke library
    int N = A.cols();									//Set N the number of the universal dimension of square matrices used
    if(B.cols() != N || A.rows() != N || B.rows() != N)					//if one of the A, B matrices is not NxN
        return false;
    int LDA = A.outerStride(); int LDB =  B.outerStride(); int LDV = Eivecs.outerStride(); int INFO = 0;

    double * alphar = Eivals.col(0).data();
    double * alphai = Eivals.col(1).data();
    double * beta = Eivals.col(2).data();
    INFO = LAPACKE_dggev(LAPACK_COL_MAJOR, 'N', 'V', N, A.data(), LDA, B.data(), LDB, alphar, alphai, beta, 0, LDV, Eivecs.data(), LDV);
    return INFO == 0;

}


bool ProblemSolver::computeL0(MatrixXd & L0)
{												//Create and compute L0
    MatrixXd identity = - MatrixXd::Identity(dimensionM, dimensionM);				//identity is -Im (Intentity matrix)
    MatrixXd zeroes = MatrixXd::Zero(dimensionM, dimensionM);					//zeroes is a matrix with only zeroes
    for (int i = 0; i < degree - 1; ++i) {                         				//for every row except the last
        L0.block(dimensionM*i,dimensionM*(i+1), dimensionM, dimensionM) = identity;		//identity is on the block matrix after the diagonal
        for (int j = 0; j < i + 1; ++j) {
            L0.block(dimensionM*i,dimensionM*j, dimensionM, dimensionM) = zeroes;		//set to zero  before the identity
        }
        for (int j = i + 2; j < degree; ++j) {
            L0.block(dimensionM*i,dimensionM*j, dimensionM, dimensionM) = zeroes;		//set to zero  after the diagonal
        }
    }
    for (int i = 0; i < degree; ++i) {								//The last row of blocks contains Mi, i= 1, 2, ...
        MatrixXd mI(dimensionM, dimensionM);
        initMI(mI, this->sylvesterPolynomial->getMd(i), dimensionM);
        L0.block(dimensionM*(degree - 1), dimensionM*i, dimensionM, dimensionM) = mI;
    }
    return true;
}

bool ProblemSolver::computeL1(MatrixXd & L1)
{												//Create and compute L1
    const int dimensionM = this->sylvesterPolynomial->getMatrixDimensions();
    MatrixXd identity = MatrixXd::Identity(dimensionM, dimensionM);				//identity is Im (Intentity matrix)
    MatrixXd zeroes = MatrixXd::Zero(dimensionM, dimensionM);					//zeroes is a matrix with only zeroes
    for (int i = 0; i < degree; ++i) {								//for every row except the last
        if(i == degree - 1)
            L1.block(dimensionM*i,dimensionM*i, dimensionM, dimensionM) = this->Md;		//if on last row diagonal is Md
        else
            L1.block(dimensionM*i,dimensionM*i, dimensionM, dimensionM) = identity;		//identity is on the block matrix diagonal
        for (int j = 0; j < i; ++j) {
            L1.block(dimensionM*i,dimensionM*j, dimensionM, dimensionM) = zeroes;		//set to zero before the identity
        }
        for (int j = i + 1; j < degree ; ++j) {
            L1.block(dimensionM*i,dimensionM*j, dimensionM, dimensionM) = zeroes;		//set to zero after the diagonal
        }
    }
    return true;
}

bool ProblemSolver::initMd(double ** tempMD, int matrixDimensions)				//Intitialize Md matrix
{
    if(tempMD == NULL)										//no matrix given
        return false;
    this->Md = MatrixXd(matrixDimensions, matrixDimensions);
    for (int i = 0; i < matrixDimensions; ++i) {
        for (int j = matrixDimensions - 1; j >= 0 ; --j) {
            this->Md(i, matrixDimensions - 1 - j) = tempMD[i][j];				//Insert given Md to Md in class by coping
        }
    }
    return true;
}

bool ProblemSolver::initMI(MatrixXd& mI, double ** tempMD, int matrixDimensions)		//Intitialize MI matrix
{
    if(tempMD == NULL)
        return false;										//no matrix given
    for (int i = 0; i < matrixDimensions; ++i) {
        for (int j = matrixDimensions - 1; j >= 0 ; --j) {
            mI(i, matrixDimensions - 1 - j) = tempMD[i][j];					//Insert given Md to Md in class by coping
        }
    }
    return true;
}

double ProblemSolver::computeStateIndicator()					//Compute k (state indicator of Md)
{
    JacobiSVD<MatrixXd> svd(Md, ComputeFullU);
    JacobiSVD<MatrixXd>::SingularValuesType svds = svd.singularValues();
//    cout<<svds<<endl;
    double max = svds(0);
    double min = svds(0);
    for (int i = 1; i < svds.rows(); ++i) {					//find max and min of singular values
        if(svds(i) > max)
            max = svds(i);
        if(svds(i) < min)
            min = svds(i);
    }
    if(min == 0) {								//is Md invertible?
        mdIsInvertible = false;
        return 1;
    }

    return max/min;								//return
}

unsigned int ProblemSolver::computeUpperBound(unsigned int power)		//compute 10^B
{
    this->upperBound = 10;
    for (int i = 0; i < power - 1; ++i) {
        this->upperBound *= 10;
    }
}
bool ProblemSolver::removeSolsWithMultiplicityStandard(Eigen::MatrixXcd & Eivecs, Eigen::MatrixXcd & Eivals, int * Multiplicity)
{						//Eliminate eigenvalues and eigenvectors with multiplicity > 1 in standard eigenproblem
    bool ret = false;
    int numRows = Eivals.rows() - 1;
    int numCols = Eivecs.cols() - 1;
    for (int j = 0; j < Eivals.rows(); ++j) {
        for (int l = 0; l < Eivals.rows(); ++l) {
            if((Eivals(j).real() == Eivals(l).real() && Eivals(j).imag() == Eivals(l).imag() && l != j)) {			//if evals(i)==eivals(j) with i =/= j
                Multiplicity[j]++;												//increase multiplicity by 1
                if( l < numRows )
                    Eivals.block(l, 0, numRows - l, Eivals.cols()) = Eivals.block(l + 1, 0, numRows - l, Eivals.cols());	//
                if( l < numCols )												//
                    Eivecs.block(0,l, Eivecs.rows(), numCols - l) = Eivecs.block(0,l+1,Eivecs.rows(),numCols-l);		//
                Eivals.conservativeResize(numRows,Eivals.cols());								//Deleting the extra eigenvalue in Eivals
                Eivecs.conservativeResize(Eivecs.rows(),numCols);								//and the extra eigenvector in Eivecs
                numRows = Eivals.rows() - 1;											//
                numCols = Eivecs.cols() - 1;											//
                ret = true;
            }
        }    
    }
    for (int j = 0; j < Eivals.rows(); ++j) {
        if(this->StandardNoSolution(Eivecs, Eivals, j, Multiplicity))								//Eigenvalue produces no solution
        {															//in standard eigenproblem
            if( j < numRows )													//
                Eivals.block(j, 0, numRows - j, Eivals.cols()) = Eivals.block(j + 1, 0, numRows - j, Eivals.cols());		//
            if( j < numCols )													//
                Eivecs.block(0,j, Eivecs.rows(), numCols - j) = Eivecs.block(0,j+1,Eivecs.rows(),numCols-j);			//Deleting the no solution eigenvalue
            Eivals.conservativeResize(numRows,Eivals.cols());									//and the eigenvector
            Eivecs.conservativeResize(Eivecs.rows(),numCols);									//
            numRows = Eivals.rows() - 1;											//
            numCols = Eivecs.cols() - 1;											//
            continue;
        }
    }
    return ret;
}

bool ProblemSolver::removeSolsWithMultiplicityGeneralized(Eigen::MatrixXd & Eivecs, Eigen::MatrixXd & Eivals, int * Multiplicity)
{						//Eliminate eigenvalues and eigenvectors with multiplicity > 1 in Generalized eigenproblem
    bool ret = false;
    int numRows = Eivals.rows() - 1;
    int numCols = Eivecs.cols() - 1;
    for (int j = 0; j < Eivals.rows(); ++j) {
        for (int l = 0; l < Eivals.rows(); ++l) {
            if(Eivals(j,0) == Eivals(l,0) && Eivals(j,2) == Eivals(l,2) && l != j ) {						//if evals(i)==eivals(j) with i =/= j
                Multiplicity[j]++;												//increase multiplicity by 1
                if( l < numRows )
                    Eivals.block(l, 0, numRows - l, Eivals.cols()) = Eivals.block(l + 1, 0, numRows - l, Eivals.cols());	//
                if( l < numCols )												//
                    Eivecs.block(0,l, Eivecs.rows(), numCols - l) = Eivecs.block(0,l+1,Eivecs.rows(),numCols-l);		//
                Eivals.conservativeResize(numRows,Eivals.cols());								//Deleting the extra eigenvalue in Eivals
                Eivecs.conservativeResize(Eivecs.rows(),numCols);								//and the extra eigenvector in Eivecs
                numRows = Eivals.rows() - 1;											//
                numCols = Eivecs.cols() - 1;											//
                ret = true;
            }
        }
    }
    for (int j = 0; j < Eivals.rows(); ++j) {
        if(this->GeneralizedNoSolution(Eivecs, Eivals, j, Multiplicity))							//Eigenvalue produces no solution
        {															//in generalized eigenproblem
            if( j < numRows )													//
                Eivals.block(j, 0, numRows - j, Eivals.cols()) = Eivals.block(j + 1, 0, numRows - j, Eivals.cols());		//
            if( j < numCols )													//
                Eivecs.block(0,j, Eivecs.rows(), numCols - j) = Eivecs.block(0,j+1,Eivecs.rows(),numCols-j);			//Deleting the no solution eigenvalue
            Eivals.conservativeResize(numRows,Eivals.cols());									//and the extra eigenvector in Eivecs
            Eivecs.conservativeResize(Eivecs.rows(),numCols);									//
            numRows = Eivals.rows() - 1;											//
            numCols = Eivecs.cols() - 1;											//
            continue;
        }
    }
    return ret;
}

bool ProblemSolver::substituteChangeOfVariable(ChangeOfVariableCoefficients * coefs)			//Change y
{
    for (int j = 0; j < this->numOfSolutions; ++j) {
        double y = this->Solutions[j]->getY();								//get x
        double x = this->Solutions[j]->getX();								//get y
        int mul = this->Solutions[j]->getMultiplicity();						//get multiplicity
        y = (coefs->t2 - y * coefs->t4)/(y * coefs->t3 - coefs->t1);					//y'=(t2-y*t4)/(y*t3-t1)
        delete this->Solutions[j];
        this->Solutions[j] = new Solution(x, y, mul);							//add x,y',multiplicity
    }
}
