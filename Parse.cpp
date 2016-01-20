#include <cstdio>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include "Parse.h"
#include <vector>
using namespace std;
using namespace Eigen;

void Parser:: readInput(char filename[], char*& polynomial1, char*& polynomial2)		//Read stdin of file
{
    if(filename == NULL)
        Parser::readFromStdin(polynomial1, polynomial2);            				//read from stdin
    else
        Parser::readFromFile(filename, polynomial1, polynomial2);       			//else read from file
}

void Parser::readPoints(MatrixXd & pointsMatrix, MatrixXd & pointsMatrix2, istream & input)
{
	string points;
	std::vector <string> points1;
	std::vector<string> points2;
	char delim = ' ';
	int num = 0;
	while(getline(cin, points, delim) != NULL)
	{
		if(!strcmp(points.c_str(),"-d"))
				break;
		points1.push_back(points);
		if(getline(cin, points, '\n') != NULL)
			points2.push_back(points);
		else
			break;
		num+=2;
	}
	num = num % 2  == 0 ? num/2 : num/2 + 1;
	pointsMatrix =  MatrixXd::Zero(num, 2);
	for(int i = 0; i < num; i++)
	{
		pointsMatrix(i, 0) = atof(points1[i].c_str());
		pointsMatrix(i, 1) = atof(points2[i].c_str());
//		cout << pointsMatrix(i, 0) << pointsMatrix(i, 1) << " " << i << endl;
	}
	points1.clear();
	points2.clear();
	num = 0;
	while(getline(cin, points, delim) != NULL)
	{
		points1.push_back(points);
		if(getline(cin, points, '\n') != NULL)
			points2.push_back(points);
		else
			break;
		num+=2;
	}
	num = num % 2  == 0 ? num/2 : num/2 + 1;
	pointsMatrix2 =  MatrixXd::Zero(num, 2);
	for(int i = 0; i < num; i++)
	{
		pointsMatrix2(i, 0) = atof(points1[i].c_str());
		pointsMatrix2(i, 1) = atof(points2[i].c_str());
//		cout << pointsMatrix2(i, 0) << pointsMatrix2(i, 1) << " " << num << " " << i << endl;
	}
}

void Parser:: readFromFile(char filename[], char*& polynomial1, char*& polynomial2)		//Read from file
{
    FILE * input;
    input=fopen(filename, "r");									//Open file
    Parser::readData(input, polynomial1, polynomial2);						//Read the data
    fclose(input);										//Close file
}

int Parser:: countLine(FILE * f)           							//count the number of chars in the line and rewind
{
    if(f==stdin)
        return 200;
    int i=0;
    char c;
    while((c=fgetc(f))!='\n')
    {
        i++;

    }
    fseek(f, -i-1, SEEK_CUR);
    return i;
}
void Parser:: readLine(FILE * f, char *s)           						//read the line
{
    int i=0;
    char c;
    while((c=fgetc(f))!='\n' && (char)c!=EOF)
    {
        s[i]=(char)c;
        i++;
    }
    s[i]='\0';
}
void Parser:: readData(FILE * input, char*& polynomial1, char*& polynomial2)			//Read data
{
    char strd1[5], strd2[5];
    polynomial1 = new char[Parser::countLine(input)+1];
    Parser:: readLine(input, polynomial1);							//read polynomial 1
    polynomial2 = new char[Parser::countLine(input)+1];
    Parser:: readLine(input, polynomial2);							//read polynomial 2
}

void Parser:: readFromStdin(char*& polynomial1, char*& polynomial2)				//Read from stdin
{
    Parser::readData(stdin, polynomial1, polynomial2);


}

