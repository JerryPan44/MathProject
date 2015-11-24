#ifndef _SOLUTION_
#define _SOLUTION_

class Solution
{
    double y;
    double x;
    int multiplicity;
public:
    Solution(double X, double Y, int Multiplicity):x(X), y(Y), multiplicity(Multiplicity)
    {

    }
    double getX()
    {
        return x;
    }
    double getY()
    {
        return y;
    }
    int getMultiplicity()
    {
        return this->multiplicity;
    }
    void PrintSolution();
};

#endif