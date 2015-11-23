#include <cstdio>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include "Parse.h"
using namespace std;

void Parser:: readInput(char filename[], char*& polynomial1, char*& polynomial2)
{
    if(filename == NULL)
        Parser::readFromStdin(polynomial1, polynomial2);            //read from stdin
    else
        Parser::readFromFile(filename, polynomial1, polynomial2);       //read from file
}

void Parser:: readFromFile(char filename[], char*& polynomial1, char*& polynomial2)
{
    FILE * input;
    input=fopen(filename, "r");
    //might change
    Parser::readData(input, polynomial1, polynomial2);
    //convert polynomial1 to matrix representation
    fclose(input);
}
int Parser:: countLine(FILE * f)            //count the number of chars in the line and rewind
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
void Parser:: readLine(FILE * f, char *s)           //read the line
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
void Parser:: readData(FILE * input, char*& polynomial1, char*& polynomial2)
{
    char strd1[5], strd2[5];
    polynomial1 = new char[Parser::countLine(input)+1];     //read polynomial 1
    Parser:: readLine(input, polynomial1);
    polynomial2 = new char[Parser::countLine(input)+1];     //read polynomial 2
    Parser:: readLine(input, polynomial2);
}

void Parser:: readFromStdin(char*& polynomial1, char*& polynomial2)
{
    //FILE * input = stdin;
    Parser::readData(stdin, polynomial1, polynomial2);


}

