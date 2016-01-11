#include "Solution.h"
#include <iostream>
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
