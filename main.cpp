#include <iostream>
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
#include <cmath>
using namespace std;
using namespace Eigen;
unsigned int getHiddenMaxDeg(SylvesterMatrix * SM, BivariatePolynomial * Bp1, BivariatePolynomial * Bp2);

void cleanResources(ProblemSolver * PS, SylvesterMatrix * , BivariatePolynomial * , BivariatePolynomial * , SylvesterPolynomial * , MyMatrix *   );
void solveProblem(char * filename, int d1, int d2,int B);
void solveGeneratedProblem(int d1, int d2, int B);
void changeOfVariable(BivariatePolynomial * Bp1, BivariatePolynomial * Bp2,SylvesterMatrix * SM, ProblemSolver * PS, int B);
bool backSubstituteSols(BivariatePolynomial * Bp1, BivariatePolynomial * Bp2, ProblemSolver * PS);

int main (int argc, char *argv[])
{
    char * filename = NULL;
    char * input;
    int d1, d2, B = 7;										//Default B=7
    bool read = false, generate = false;
    if(argc > 10 || argc < 7)
    {
        perror("Wrong Number of agruments");							//Failed input arguments
        return 0;
    }
    for (int i = 0; i < argc; ++i) {
        if(!strcmp(argv[i], "-read"))
            read = true;									//Read from input
        if(read) {
            if (!strcmp(argv[i], "-i")) {
                filename = argv[i + 1];								//filename after -i
                if (access(filename, F_OK) == -1)        					//if file doesnt exist
                {
                    perror("Error : File does not exist");
                    return 0;
                }
            }
        }
        if(!strcmp(argv[i], "-d1"))
            d1 = atoi(argv[i+1]);								//d1 from input after -d1
        if(!strcmp(argv[i], "-d2"))
            d2 = atoi(argv[i+1]);								//d2 from input after -d2
        if(!strcmp(argv[i], "-solve"))
            B = atoi(argv[i+1]);								//B from input after -solve
        if(!strcmp(argv[i], "-generate"))
            generate = true;									//generate from input
    }

    if(generate)
        solveGeneratedProblem(d1 , d2, B);							//Solve the problem with generate
    else
        solveProblem(filename, d1, d2, B);							//Else solve with given polynomials

    return 0;											//End of main

}


void solveGeneratedProblem(int d1, int d2, int B)						//Generate and Solve
{
    BivariatePolynomial * Bp1 = new BivariatePolynomial(d1);
    BivariatePolynomial * Bp2 = new BivariatePolynomial(d2);
    SylvesterMatrix * SM = new SylvesterMatrix(Bp1, Bp2);
    int SpDeg = getHiddenMaxDeg(SM, Bp1, Bp2);							//Which is the hidden variable?
    SylvesterPolynomial * SP = new SylvesterPolynomial(SpDeg, SM->getRowDimension());
    SP->SMatrixToSPolynomial(SM);
    MyMatrix * m = new MyMatrix(SM->getRowDimension(), 1);      				//generate random matrix of 1 column (vector v)
    ProblemSolver * PS = new ProblemSolver(SP, B, SP->getDegree());				//create a problem solver class
    PS->Solve();										//solve
    backSubstituteSols(Bp1, Bp2, PS);								//replace sols and check if they are indeed sols
    SM->changeOfVariable();									//change of variable logic
    changeOfVariable(Bp1, Bp2, SM, PS, B);
    cleanResources(PS, SM, Bp1, Bp2, SP, m);        						//clean resources
}


void solveProblem(char * filename, int d1, int d2,int B)
{
    char * polynomial1, * polynomial2;
    Parser::readInput(filename, polynomial1, polynomial2);					//read d1, d2, p1, p2
    fprintf(stdout,"Equations \n------- \n%s\n%s\n ", polynomial1, polynomial2);
    int systemDegree=0;
    BivariatePolynomial * Bp1 = new BivariatePolynomial(polynomial1, d1);			//initialize Bp1, Bp2
    BivariatePolynomial * Bp2 = new BivariatePolynomial(polynomial2, d2);
    SylvesterMatrix * SM = new SylvesterMatrix(Bp1, Bp2);					//construct the sylvester matrix from Bp1 and Bp2
    int SpDeg = getHiddenMaxDeg(SM, Bp1, Bp2);
    SylvesterPolynomial * SP = new SylvesterPolynomial(SM->getHiddenDeg(), SM->getRowDimension());	//sylvester polynomial of sylvester matrix hidden variable degree
    SP->SMatrixToSPolynomial(SM);								//convert Sylvester Matrix to sylvester polynomial
    MyMatrix * m = new MyMatrix(SM->getRowDimension(), 1);					//random matrix of 1 column (the vector v)
/*  m->generate();
    cout<<"multiplying Sylvester Matrix with"<<endl;
    m->Print();
    SylvesterMatrix * result;
    SM->multiply(m, result);                    						//multiply the sylvester matrix with the random vector v
    result->Print();                           							//print the result
    //Project 2 code */
    ProblemSolver * PS = new ProblemSolver(SP, B, SP->getDegree());
    PS->Solve();
    backSubstituteSols(Bp1, Bp2, PS);
    changeOfVariable(Bp1, Bp2, SM, PS, B);
//    delete result;
    cleanResources(PS, SM, Bp1, Bp2, SP, m);							//clean resources
}

bool backSubstituteSols(BivariatePolynomial * Bp1, BivariatePolynomial * Bp2, ProblemSolver * PS)//Compute polynomials for computed solutions
{
    int numSols = PS->getNumOfSols();
    for (int i = 0; i < numSols; ++i) {								//For every solution
        Solution * sol = PS->getSolution(i);
        if(sol->getMultiplicity() == 1) {                           //if multiplicity of solution is 1 back substitute the sols and check if they nullify the polynomials
            double res1 = Bp1->backSubstitute(sol->getX(), sol->getY());			//Compute the value of the polynomial Bp1
            double res2 = Bp2->backSubstitute(sol->getX(), sol->getY());			//Compute the value of the polynomial Bp2
            if(abs(res1) < 0.000001 && abs(res2) < 0.000001)					//Are the values over 10^-6
                cout<<endl<<"SOLUTION : y = "<<sol->getY()<<" x = "<<sol->getX()<<" ACCEPTED"<<endl;
            else
                cout<<endl<<"SOLUTION : y = "<<sol->getY()<<" x = "<<sol->getX()<<" REJECTED : "<<res1<<" , "<<res2<<endl;
        }
    }
}
void changeOfVariable(BivariatePolynomial * Bp1, BivariatePolynomial * Bp2, SylvesterMatrix * SM, ProblemSolver * PS, int B)
{												//Change of variable as asked in 4)
    SylvesterPolynomial * SP;
    ProblemSolver * PSNew;
    for (int i = 0; i < 3; ++i) {                                       //try to change variable with random t1,t2,t3,t4 up to 3 times until you find a better k
        SylvesterMatrix * SmTemp = new SylvesterMatrix(*SM);
        SP = new SylvesterPolynomial(SmTemp->getHiddenDeg(), SmTemp->getRowDimension());
        SmTemp->changeOfVariable();								//Pick 4 random t(i) in Sylvester marix
        SP->SMatrixToSPolynomial(SmTemp);                           				//Convert Sylvester Matrix to sylvester polynomial
//        SP->Print();
        PSNew = new ProblemSolver(SP, B , SP->getDegree());					//Compute Md', k'(state indicator)
//        cout<<"OLD K : "<<PS->getStateIndicator()<<" NEW K : "<<PSNew->getStateIndicator()<<endl<<endl;
        if(PS->getStateIndicator() > PSNew->getStateIndicator()
           || (PS->getMdIsInvertible() == false && PSNew->getMdIsInvertible() == true))				//Check for variable change (k(Md) > k'(Md'))
        {
            cout<<"************SOLVING WITH CHANGE OF VARIABLE**************"<<endl;
            PSNew->Solve();									//Solve the new eigenproblem
            ChangeOfVariableCoefficients * coefs = SmTemp->getCoefs();
            PSNew->substituteChangeOfVariable(coefs);						//compute new y
            backSubstituteSols(Bp1, Bp2, PSNew);						//Compute polynomials for computed solutions
            delete PSNew;
            delete SP;
            return;
        }
        delete PSNew;
        delete SP;
    }
}
void cleanResources(ProblemSolver * PS, SylvesterMatrix * SM, BivariatePolynomial * Bp1, BivariatePolynomial * Bp2, SylvesterPolynomial * SP, MyMatrix * m)
{												//Delete created resources
    delete PS;
    delete SM;
    delete Bp1;
    delete Bp2;
    delete SP;
    delete m;
}



unsigned int getHiddenMaxDeg(SylvesterMatrix * SM, BivariatePolynomial * Bp1, BivariatePolynomial * Bp2)
{
    if(SM->getHiddenVariable() == 'y')
    {
        return Bp1->getdegy() > Bp2->getdegy()? Bp1->getdegy() : Bp2->getdegy();
    }
    if(SM->getHiddenVariable() == 'x')
    {
        return Bp1->getdegx() > Bp2->getdegx()? Bp1->getdegx() : Bp2->getdegx();
    }

}
