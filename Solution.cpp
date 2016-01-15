#include "Solution.h"
#include <iostream>
#include <cmath>
#include "ErrorMargin.h"
using namespace std;
void Solution::PrintSolution()							//Print y, x, multiplicity
{
    cout<<"y = "<<y<<",";
    if(this->multiplicity > 1)
    {
    	cout<<" multiplicity = "<<multiplicity;
    	for(int i = 0 ; i < this->lastSet + 1; i++)
    	{
    		cout<<" x = "<<x[i];
    	}
    	cout<<endl;
    }
    else
        cout<<" x = "<<x[0]<<endl;
}

void Solution::PrintSolution(ofstream & file)							//Print y, x, multiplicity
{
    if(this->multiplicity > 1)
    {
        for(int i = 0 ; i < this->lastSet + 1; i++)
        {
            file << this->x[i] << " ";
            file << this->y << " 0" << endl;
        }
    }
    else
    {
        file << this->x[0] << " ";
        file << this->y << " 0" << endl;
    }

}

bool Solution::compareSolutions(double x, double y)
{
    if(this->multiplicity > 1)
    {
        for(int i = 0 ; i < this->lastSet + 1; i++)
        {
            if(abs(x - this->x[i]) < ERROR_MARGIN && abs(y - this->y) < ERROR_MARGIN)
                return true;
        }
    }
    else
    {
        if(abs(x - this->x[0]) < ERROR_MARGIN && abs(y - this->y) < ERROR_MARGIN)
            return true;
    }
    return false;
}
