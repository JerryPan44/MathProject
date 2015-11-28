#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <time.h>
#include <cstring>
#include "BivariatePolynomial.h"
using namespace std;

void BivariatePolynomial::InitDegs()							//find polynomial y, x degrees
{
    this->degy = 0;
    bool foundx=false,foundy=false;
    this->degx = 0;
    for (int i = matrixDimension -1; i >=0 ; --i) {
        for(int j=matrixDimension-1;j>=0 ; j--)
            if(this->MatrixRepresentation[i][j]!=0)					//find the biggest power of y
            {
                this->degy = i;
                foundx = true;
                break;
            }
        if(foundx)
            break;
    }
    for (int i = matrixDimension-1; i >=0 ; --i) {					//find the biggest power of x
        for(int j=matrixDimension-1;j>=0 ; j--)
            if(this->MatrixRepresentation[j][i]!=0)
            {
                foundy = true;
                this->degx = i;
                break;
            }
        if(foundy)
            break;
    }

}

void BivariatePolynomial::Print()							//Print polynomial
{
    cout<<"DEGY : "<<this->degy<<endl<<"DEGX : "<<this->degx<<endl;			//print degree of y, degree of x
    for(int i=matrixDimension-1;i >= 0;i--)
    {
        for (int j = matrixDimension - 1; j >= 0; --j)
        {
            cout << this->MatrixRepresentation[i][j]<<" ";
        }
        cout<<endl;
    }
}

BivariatePolynomial::BivariatePolynomial(char polstr[], int sD)
{
    matrixDimension = sD + 1;
    this->MatrixRepresentation = new int*[matrixDimension];
    for (int i = 0; i < matrixDimension; ++i) {
        this->MatrixRepresentation[i] = new int[matrixDimension];
        for (int j = 0; j < matrixDimension; ++j) {
            this->MatrixRepresentation[i][j] = 0;					//initialize matrix elements (== 0)
        }
    }
    Parse(polstr);									//parse the polynomial
    InitDegs();										//find the degrees
}
BivariatePolynomial::BivariatePolynomial(int sD)					//Create random polynomial
{

    srand (time(NULL)+rand()%10000);
    matrixDimension = sD + 1;
    this->MatrixRepresentation = new int*[matrixDimension];
    for (int i = 0; i < matrixDimension; ++i) {
        this->MatrixRepresentation[i] = new int[matrixDimension];
        for (int j = 0; j < matrixDimension; ++j) {
            this->MatrixRepresentation[i][j] = (rand()%100) - 50;			//Matrix element range (-50,50)
        }
    }
    InitDegs();										//Find the degrees
}


void BivariatePolynomial::Parse(char polstr[])
{
    if(this->MatrixRepresentation == NULL)
    {
        perror("Matrix is null");							//Empty Matrix
        return;
    }
    int col=0,row=0,coef=1,j=0;
    bool isPositive = true;

    if(polstr[0] == '-') {                              				//zero power coeficient case
        isPositive = false;
        j++;
    }
    if(isPositive)
        this->MatrixRepresentation[0][0] = getCoefFromPolynomial(j, polstr);
    else
        this->MatrixRepresentation[0][0] = -getCoefFromPolynomial(j, polstr);

    int choice = 0;
    int i = j;
    while(i <= strlen(polstr)){
        col = 0;
        row = 0;
        if (polstr[i] == '+' || polstr[i] == '-') {
            isPositive = polstr[i] == '+' ? true : false;           			//find if coef is positive or negative
            i++;
            coef = getCoefFromPolynomial(i, polstr);    				//reads the coeficient
            choice = 0;
            ParsePolynomialBody(i, polstr, 'x', 'y', col, row, choice, coef);           //parse either with *x first
            ParsePolynomialBody(i, polstr, 'y', 'x', row, col, choice, coef);           //or with *y first
            if (isPositive)
                this->MatrixRepresentation[row][col] = coef;
            else
                this->MatrixRepresentation[row][col] = -coef;
        }
        else
        {
            i++;
        }

    }

}

void BivariatePolynomial::ParsePolynomialBody(int & i, char * polstr, char first, char second, int & firstCoord, int & secondCoord, int & choice, int coef)
											//parsing method
{
    if(choice == 1)									//if you chose another method before that return
        return;
    if(polstr[i + 1] == second || (polstr[i + 2] == second && coef == 1))		//case we arent in the correct call
        return;

    if (polstr[i + 1] == '*' || coef == 1) {
        if(coef == 1 && polstr[i+1] != '*')						//case we have x*y and not 1*x*y
            i--;
        if (polstr[i + 2] == first) {
            choice = 1;
            if (polstr[i + 3] == '^') {
                if(!consume(polstr, i, 4))
                    return;								//no of chars consumed
                firstCoord = getCoefFromPolynomial(i, polstr);
                if (polstr[i + 1] == '*') {
                    if (polstr[i + 2] == second) {
                        if (polstr[i + 3] == '^') {
                            if(!consume(polstr, i, 4))
                                return;
                            secondCoord = getCoefFromPolynomial(i, polstr);
                        }
                        else
                        {
                            secondCoord = 1;
                            if(!consume(polstr, i, 3))
                                return;
                        }
                    }
                }
            }
            else {
                if(!consume(polstr, i, 3))
                    return;
                firstCoord = 1;
                if (polstr[i] == '*') {
                    if (polstr[i + 1] == second) {
                        if (polstr[i + 2] == '^') {
                            if(!consume(polstr, i, 3))
                                return;
                            secondCoord = getCoefFromPolynomial(i, polstr);		//no of chars consumed
                        }
                        else
                        {
                            secondCoord = 1;
                            if(!consume(polstr, i, 2))
                                return;
                        }
                    }
                }
            }
        }
    }
}

bool BivariatePolynomial::consume(char * str, int & pos, int numChars)
{
    for (int i = 0; i < numChars; ++i) {
        if(str[pos] == '\0')
            return false;
        pos++;
    }
    return true;
}

int BivariatePolynomial::getCoefFromPolynomial(int & pos, char * polstr)
{
    int count = 0;
    if(!isdigit(polstr[pos]))
    {
        pos += count - 1;
        return 1;									//if no digit return 1
    }

    char * strCoef;
    while (isdigit(polstr[pos + count])) {
        count++;
    }
    strCoef = new char[count+1];
    strCoef[count] = '\0';
    count = 0;
    while (isdigit(polstr[pos + count])) {						//else read the digits one by one and return them as an int
        strCoef[count] = polstr[pos + count];
        count++;
    }
    pos +=count-1;
    int k = atoi(strCoef);
    delete []strCoef;
    return k;
}
BivariatePolynomial::~BivariatePolynomial()
{
    delete []this->MatrixRepresentation;
}

double BivariatePolynomial::exp(double num, double power)
{
    if(power == 0)
        return 1;
    double first = num;
    for (int i = 0; i < power - 1; ++i) {
        num *= first;
    }
    return num;
}

double BivariatePolynomial::backSubstitute(double x, double y)
{
    double sum = 0;
    for (int i = 0; i < this->degy + 1; ++i) {
        for (int j = 0; j < this->degx + 1; ++j) {
            if(this->MatrixRepresentation[i][j] != 0)
                sum += this->MatrixRepresentation[i][j]*(exp(y, i) * exp(x, j));
        }
    }
    return sum;
}
