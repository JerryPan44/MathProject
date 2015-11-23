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
    void PrintSolution();
};

#endif