#include "Solution.h"
#include <iostream>
using namespace std;
void Solution::PrintSolution()							//Print y, x, multiplicity
{
    cout<<"y = "<<y<<",";
    if(this->multiplicity > 1)
        cout<<" multiplicity = "<<multiplicity<<endl;
    else
        cout<<" x = "<<x<<endl;
}
