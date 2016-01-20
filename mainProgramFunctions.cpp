#include "mainProgramFunctions.h"

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
//    backSubstituteSols(Bp1, Bp2, PS);								//replace sols and check if they are indeed sols
    SM->changeOfVariable();									//change of variable logic
    changeOfVariable(Bp1, Bp2, SM, PS, B);
    cleanResources(PS, SM, Bp1, Bp2, SP, m);        						//clean resources
}


void solveProblem(char * filename, int d1, int d2,int B, bool fromPoints)
{
    char * polynomial1, * polynomial2;
    BivariatePolynomial * Bp1;
    BivariatePolynomial * Bp2;
    if(fromPoints == true)
    {
        MatrixXd pointsMatrix, pointsMatrix2;
        Parser::readPoints(pointsMatrix, pointsMatrix2, cin);
        cout<<"Points 1 :"<<endl;
        cout<<pointsMatrix<<endl;
        cout<<"Points 2 :"<<endl;
        cout<<pointsMatrix2<<endl;
        ofstream equationsTxt;
        equationsTxt.open("InterpolationEquations.txt");
        Interpolation interpolation(d1, pointsMatrix), interpolation2(d2 , pointsMatrix2);
        Bp1 = interpolation.find(equationsTxt);
        Bp2 = interpolation2.find(equationsTxt);
        if(Bp1 == NULL || Bp2 == NULL) {
            cout << "Error in interpolation" << endl;
            exit(0);
        }
        equationsTxt.close();
        //interpolation logic
    }
    else
    {
        Parser::readInput(filename, polynomial1, polynomial2);					//read d1, d2, p1, p2
        fprintf(stdout,"Equations \n------- \n%s\n%s\n ", polynomial1, polynomial2);
        int systemDegree=0;
        Bp1 = new BivariatePolynomial(polynomial1, d1);			//initialize Bp1, Bp2
        Bp2 = new BivariatePolynomial(polynomial2, d2);
    }
    SylvesterMatrix * SM = new SylvesterMatrix(Bp1, Bp2);					//construct the sylvester matrix from Bp1 and Bp2
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
    findFullSolutions(Bp1, Bp2, PS);
    changeOfVariable(Bp1, Bp2, SM, PS, B);
//    delete result;
    cleanResources(PS, SM, Bp1, Bp2, SP, m);							//clean resources
}

void findFullSolutions(BivariatePolynomial * Bp1, BivariatePolynomial * Bp2, ProblemSolver * PS)
{
    int numSols = PS->getNumOfSols();
    for (int i = 0; i < numSols; ++i) {								//For every solution
        Solution * sol = PS->getSolution(i);
        if(sol->getMultiplicity() > 1) {                           //if multiplicity of solution is 1 back substitute the sols and check if they nullify the polynomials
            Polynomial * result1 = Bp1->backSubstitute(sol->getY(), PS->getHiddenVar());
            Polynomial * result2 = Bp2->backSubstitute(sol->getY(), PS->getHiddenVar());
            int numRoots1, numRoots2;
            double * roots1 = result1->computeAndGetRoots(numRoots1);
            double * roots2 = result2->computeAndGetRoots(numRoots2);
            int index = 0;
//            cout<<"roots 1"<<endl;
//            for(int i = 0; i < numRoots1; i++)
//                cout<< roots1[i]<<endl;
//            cout<<"roots 2"<<endl;
//            for(int i = 0; i < numRoots2; i++)
//                cout<< roots2[i]<<endl;
            for(int i = 0; i < numRoots1; i++)
            {
                for(int j = 0; j < numRoots2; j++)
                {
                    if(abs(roots1[i] - roots2[j]) < ERROR_MARGIN)
                    {
                        sol->setX(roots1[i], index);
                        index++;
                    }
                }
            }
        }
//        sol->PrintSolution();
    }
}


int getMax(int m1, int m2)
{
    if(m1>m2)
        return m1;
    return m2;
}

bool backSubstituteSols(BivariatePolynomial * Bp1, BivariatePolynomial * Bp2, ProblemSolver * PS, ofstream & solutionsTxt)//Compute polynomials for computed solutions
{
    int numSols = PS->getNumOfSols();
    cout<<"Regular Roots"<<endl;
    cout<<"--------------"<<endl;
    for (int i = 0; i < numSols; ++i) {								//For every solution
        Solution * sol = PS->getSolution(i);
        if(sol->getMultiplicity() == 1) {                           //if multiplicity of solution is 1 back substitute the sols and check if they nullify the polynomials
            double res1 = Bp1->backSubstituteXandY(sol->getX(), sol->getY());			//Compute the value of the polynomial Bp1
            double res2 = Bp2->backSubstituteXandY(sol->getX(), sol->getY());			//Compute the value of the polynomial Bp2
            if(abs(res1) < ERROR_MARGIN && abs(res2) < ERROR_MARGIN)					//Are the values over 10^-6
            {
                solutionsTxt << sol->getX() << " " << sol->getY() << " 0"<< endl;
                cout<<endl<<"SOLUTION : y = "<<sol->getY()<<" x = "<<sol->getX()<<" ACCEPTED"<<endl;
            }
            else
                cout<<endl<<"SOLUTION : y = "<<sol->getY()<<" x = "<<sol->getX()<<" REJECTED : "<<res1<<" , "<<res2<<endl;
        }
        else
        if(sol->getMultiplicity() > 1)
        {
            for (int j = 0; j < sol->getLastSet() + 1; ++j) {
                double res1 = Bp1->backSubstituteXandY(sol->getXat(j), sol->getY());			//Compute the value of the polynomial Bp1
                double res2 = Bp2->backSubstituteXandY(sol->getXat(j), sol->getY());			//Compute the value of the polynomial Bp2
                if(abs(res1) < ERROR_MARGIN && abs(res2) < ERROR_MARGIN)					//Are the values over 10^-6
                {
                    solutionsTxt << sol->getXat(j) << " " << sol->getY() << " 0"<< endl;
                    cout<<endl<<"SOLUTION : y = "<<sol->getY()<<" x = "<<sol->getXat(j)<<" ACCEPTED"<<endl;
                }
                else
                    cout<<endl<<"SOLUTION : y = "<<sol->getY()<<" x = "<<sol->getXat(j)<<" REJECTED : "<<res1<<" , "<<res2<<endl;
            }
        }
    }
}

bool backSubstituteSols(BivariatePolynomial * Bp1, BivariatePolynomial * Bp2, ProblemSolver * PS, ProblemSolver * PSNew)//Compute polynomials for computed solutions
{
    ofstream solutionsTxt;
    solutionsTxt.open("solutions.txt");
    int numSols = PSNew->getNumOfSols();
    int oldNumSols = PS->getNumOfSols();
    backSubstituteSols(Bp1, Bp2, PS, solutionsTxt);
    cout<<endl<<"Roots with change of variable"<<endl;
    cout<<"--------------"<<endl;
    for (int i = 0; i < numSols; ++i) {								//For every solution
        Solution * sol = PSNew->getSolution(i);
        if(sol->getMultiplicity() == 1) {                           //if multiplicity of solution is 1 back substitute the sols and check if they nullify the polynomials
            double res1 = Bp1->backSubstituteXandY(sol->getX(), sol->getY());			//Compute the value of the polynomial Bp1
            double res2 = Bp2->backSubstituteXandY(sol->getX(), sol->getY());			//Compute the value of the polynomial Bp2
            if(abs(res1) < ERROR_MARGIN && abs(res2) < ERROR_MARGIN)					//Are the values over 10^-6
            {
                bool flag = false;
                for (int j = 0; j < oldNumSols; ++j){
                    Solution * oldSol = PS->getSolution(j);
                    if(oldSol->compareSolutions(sol->getX(), sol->getY()))
                    {
                        flag = true;
                        break;
                    }
                }
                if(!flag)
                {
                    solutionsTxt << sol->getX() << " " << sol->getY() << " 0"<< endl;
                }
                cout<<endl<<"SOLUTION : y = "<<sol->getY()<<" x = "<<sol->getX()<<" ACCEPTED"<<endl;
            }
            else
                cout<<endl<<"SOLUTION : y = "<<sol->getY()<<" x = "<<sol->getX()<<" REJECTED : "<<res1<<" , "<<res2<<endl;
        }
        else
        if(sol->getMultiplicity() > 1)
        {
            for (int j = 0; j < sol->getLastSet() + 1; ++j) {
                double res1 = Bp1->backSubstituteXandY(sol->getXat(j), sol->getY());			//Compute the value of the polynomial Bp1
                double res2 = Bp2->backSubstituteXandY(sol->getXat(j), sol->getY());			//Compute the value of the polynomial Bp2
                if(abs(res1) < ERROR_MARGIN && abs(res2) < ERROR_MARGIN)					//Are the values over 10^-6
                {
                    bool flag = false;
                    int pos;
                    for (int i = 0; i < oldNumSols; ++i){
                        Solution * oldSol = PS->getSolution(i);
                        if(oldSol->compareSolutions(sol->getXat(j), sol->getY()))
                        {
                            pos = j;
                            flag = true;
                            break;
                        }
                    }
                    if(!flag)
                    {
                        solutionsTxt << " " << sol->getXat(pos) << sol->getY() << " 0"<< endl;
                    }
                    cout<<endl<<"SOLUTION : y = "<<sol->getY()<<" x = "<<sol->getXat(j)<<" ACCEPTED"<<endl;
                }
                else
                    cout<<endl<<"SOLUTION : y = "<<sol->getY()<<" x = "<<sol->getXat(j)<<" REJECTED : "<<res1<<" , "<<res2<<endl;
            }
        }
    }
    solutionsTxt.close();
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
            backSubstituteSols(Bp1, Bp2, PS, PSNew);						//Compute polynomials for computed solutions
            delete PSNew;
            delete SP;
            return;
        }
        delete PSNew;
        delete SP;
    }
    ofstream solutionsTxt;
    solutionsTxt.open("solutions.txt");
    backSubstituteSols(Bp1, Bp2, PS, solutionsTxt);
    solutionsTxt.close();
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