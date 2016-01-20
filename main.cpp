#include <iostream>
#include "mainProgramFunctions.h"
using namespace std;
using namespace Eigen;

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
    bool fromPoints = false;
    for (int i = 0; i < argc; ++i) {
    	if(!strcmp(argv[i],"-points"))
    	{
    		fromPoints = true;
    		if(read == true)
    		{
    			perror("error cant read from file and from points");
    			return 0;
    		}
    	}
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
        solveProblem(filename, d1, d2, B, fromPoints);							//Else solve with given polynomials

    return 0;											//End of main

}


