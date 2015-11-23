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
using namespace std;
using namespace Eigen;
unsigned int getHiddenMaxDeg(SylvesterMatrix * SM, BivariatePolynomial * Bp1, BivariatePolynomial * Bp2);

void cleanResources(ProblemSolver * PS, SylvesterMatrix * , BivariatePolynomial * , BivariatePolynomial * , SylvesterPolynomial * , MyMatrix *   );
void solveProblem(char * filename, int d1, int d2,int B);
void solveGeneratedProblem(int d1, int d2, int B);
void changeOfVariable(SylvesterMatrix * SM, ProblemSolver * PS, int B);
int main (int argc, char *argv[])
{
    char * filename = NULL;
//    if(argc != 2)
//        return 1;
//    if(argc == 2)
//    {
//        filename = new char[strlen(argv[1]) + 1];
//        strcpy(filename, argv[1]);
//        if( access( filename, F_OK ) == -1 )        //if file doesnt exist
//        {
//            perror("Error in arguments");
//            return 0;
//        }
//    }
    char * input;
    int d1, d2, B = 7;
    bool read = false, generate = false;
    if(argc > 10 || argc < 7)
    {
        perror("Wrong Number of agruments");
        return 0;
    }
    for (int i = 0; i < argc; ++i) {
        if(!strcmp(argv[i], "-read"))
            read = true;
        if(read) {
            if (!strcmp(argv[i], "-i")) {
                filename = argv[i + 1];
                if (access(filename, F_OK) == -1)        //if file doesnt exist
                {
                    perror("Error : File does not exist");
                    return 0;
                }
            }
        }
        if(!strcmp(argv[i], "-d1"))
            d1 = atoi(argv[i+1]);
        if(!strcmp(argv[i], "-d2"))
            d2 = atoi(argv[i+1]);
        if(!strcmp(argv[i], "-solve"))
            B = atoi(argv[i+1]);
        if(!strcmp(argv[i], "-generate"))
            generate = true;
    }

    if(generate)
        solveGeneratedProblem(d1 , d2, B);
    else
        solveProblem(filename, d1, d2, B);

    return 0;

}


void solveGeneratedProblem(int d1, int d2, int B)
{
    BivariatePolynomial * Bp1 = new BivariatePolynomial(d1);
    //Bp1->Print();
    BivariatePolynomial * Bp2 = new BivariatePolynomial(d2);
    //Bp2->Print();
    SylvesterMatrix * SM = new SylvesterMatrix(Bp1, Bp2);
    //SM->Print();
    int SpDeg = getHiddenMaxDeg(SM, Bp1, Bp2);
    SylvesterPolynomial * SP = new SylvesterPolynomial(SpDeg, SM->getRowDimension());
    SP->SMatrixToSPolynomial(SM);
    //SP->Print();
    MyMatrix * m = new MyMatrix(SM->getRowDimension(), 1);      //random matrix of 1 column (the vector v)
    m->generate();
    SylvesterMatrix * result;
    SM->multiply(m, result);                    //multiply the sylvester matrix with the random vector v
    ProblemSolver * PS = new ProblemSolver(SP, B, SP->getDegree());
    PS->Solve();
    SM->changeOfVariable();
    changeOfVariable(SM, PS, B);
    delete result;
    cleanResources(PS, SM, Bp1, Bp2, SP, m);        //clean resources
}


void solveProblem(char * filename, int d1, int d2,int B)
{
    char * polynomial1, * polynomial2;
    Parser::readInput(filename, polynomial1, polynomial2);          //read d1 d2 p1 p2
    fprintf(stdout,"Equations \n------- \n%s\n%s\n ", polynomial1, polynomial2);
    int systemDegree=0;
    BivariatePolynomial * Bp1 = new BivariatePolynomial(polynomial1, d1);            //initialize Bp1, Bp2
    Bp1->Print();
    BivariatePolynomial * Bp2 = new BivariatePolynomial(polynomial2, d2);
    Bp2->Print();
    SylvesterMatrix * SM = new SylvesterMatrix(Bp1, Bp2);                       //construct the sylvester matrix from Bp1 and Bp2
    //SM->Print();
    int SpDeg = getHiddenMaxDeg(SM, Bp1, Bp2);
    SylvesterPolynomial * SP = new SylvesterPolynomial(SM->getHiddenDeg(), SM->getRowDimension());               //sylvester polynomial of sylvester matrix hidden variable degree
    SP->SMatrixToSPolynomial(SM);                           //convert Sylvester Matrix to sylvester polynomial
    SP->Print();
    MyMatrix * m = new MyMatrix(SM->getRowDimension(), 1);      //random matrix of 1 column (the vector v)
    m->generate();
    //cout<<"multiplying Sylvester Matrix with"<<endl;
    // m->Print();
    SylvesterMatrix * result;
    SM->multiply(m, result);                    //multiply the sylvester matrix with the random vector v
    //result->Print();                            //print the result
    //Project 2 code
    ProblemSolver * PS = new ProblemSolver(SP, B, SP->getDegree());
    PS->Solve();
    changeOfVariable(SM, PS, B);
    ProblemSolver * PsChanged = new ProblemSolver(SP, B, SP->getDegree());
    delete result;
    cleanResources(PS, SM, Bp1, Bp2, SP, m);        //clean resources
}
void changeOfVariable(SylvesterMatrix * SM, ProblemSolver * PS, int B)
{
    SylvesterPolynomial * SP;
    ProblemSolver * PSNew;
    for (int i = 0; i < 3; ++i) {
        SylvesterMatrix * SmTemp = new SylvesterMatrix(*SM);
        SP = new SylvesterPolynomial(SmTemp->getHiddenDeg(), SmTemp->getRowDimension());
        SmTemp->changeOfVariable();
        SP->SMatrixToSPolynomial(SmTemp);                           //convert Sylvester Matrix to sylvester polynomial
        SP->Print();
        PSNew = new ProblemSolver(SP, B , SP->getDegree());
        cout<<"OLD K : "<<PS->getStateIndicator()<<" NEW K : "<<PSNew->getStateIndicator()<<endl<<endl;
        if(PS->getStateIndicator() > PSNew->getStateIndicator())
        {
            cout<<"************SOLVING**************"<<endl;
            PSNew->Solve();
            delete PSNew;
            delete SP;
            return;
        }
        delete PSNew;
        delete SP;
    }
}
void cleanResources(ProblemSolver * PS, SylvesterMatrix * SM, BivariatePolynomial * Bp1, BivariatePolynomial * Bp2, SylvesterPolynomial * SP, MyMatrix * m)
{
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