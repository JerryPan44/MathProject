#ifndef _SOLUTION_
#define _SOLUTION_

#include <fstream>
#include <iostream>
#include <assert.h>

class Solution
{
    double y;
    double * x;
    int multiplicity;
    int lastSet;
public:
    Solution(double X, double Y, int Multiplicity): y(Y), multiplicity(Multiplicity)
    {
    	if(multiplicity >= 1)
    		x = new double[multiplicity];
    	x[0] = X;
    	lastSet = -1;
    }

    void setX(double  Xs, int index)
    {
    	assert(index < multiplicity);
    	this->x[index] = Xs;
    	lastSet = index;
    }

    double getY()
    {
        return y;
    }
    double getX()
    {
        return x[0];
    }
    double getXat(int i)
    {
        assert(i < (this->lastSet + 1));
        return x[i];
    }

    int getMultiplicity()
    {
        return this->multiplicity;
    }
    int getLastSet()
    {
        return this->lastSet;
    }
    void PrintSolution();
    void PrintSolution(std::ofstream &);

    bool compareSolutions(double, double);

};

#endif
