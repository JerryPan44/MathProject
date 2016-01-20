#include "../../mainProgramFunctions.h"
#include <iostream>
#include <vector>
#include <fstream>
using namespace Eigen;
using namespace std;
class MainProgramTests : public ::testing::Test
{
protected:
    virtual void SetUp() {
    }
};

void readSolutions(MatrixXd & solutions)
{
    ifstream inputFile;
    inputFile.open("solutions.txt");
    string points;
    std::vector <string> points1;
    std::vector<string> points2;
    char delim = ' ';
    int num = 0;
    while(getline(inputFile, points, delim) != NULL)
    {
        points1.push_back(points);
        if(getline(inputFile, points, delim) != NULL)
            points2.push_back(points);
        else
            break;
        if(getline(inputFile, points, '\n') == NULL)
            break;
        num+=2;
    }
    num = num % 2  == 0 ? num/2 : num/2 + 1;
    solutions =  MatrixXd::Zero(num, 2);
    for(int i = 0; i < num; i++)
    {
        solutions(i, 0) = atof(points1[i].c_str());
        solutions(i, 1) = atof(points2[i].c_str());
//		cout << pointsMatrix(i, 0) << pointsMatrix(i, 1) << " " << i << endl;
    }
}

TEST_F(MainProgramTests, TestExample1)
{
char example1[50] = "../../Examples/FromFileExamples/Example1";
std::cout.setstate(std::ios_base::failbit);
solveProblem(example1, 2, 2, 7, false);
std::cout.clear();
MatrixXd solutions;
readSolutions(solutions);
cout<<endl<<"Example1 solutions"<<endl;
cout<<"---------------------"<<endl;
cout << solutions<< endl;
}

TEST_F(MainProgramTests, TestExample2)
{
char example[50] = "../../Examples/FromFileExamples/Example2";
std::cout.setstate(std::ios_base::failbit);
solveProblem(example, 3, 3, 7, false);
std::cout.clear();
MatrixXd solutions;
readSolutions(solutions);
cout<<endl<<"Example2 solutions"<<endl;
cout<<"---------------------"<<endl;
cout << solutions<< endl;
}
TEST_F(MainProgramTests, TestExample3)
{
char example[50] = "../../Examples/FromFileExamples/Example3";
std::cout.setstate(std::ios_base::failbit);
solveProblem(example, 2, 3, 7, false);
std::cout.clear();
MatrixXd solutions;
readSolutions(solutions);
cout<<endl<<"Example 3 solutions"<<endl;
cout<<"---------------------"<<endl;
cout << solutions<< endl;
}

TEST_F(MainProgramTests, TestExample4)
{
char example[50] = "../../Examples/FromFileExamples/Example4";
std::cout.setstate(std::ios_base::failbit);
solveProblem(example, 2, 3, 7, false);
std::cout.clear();
MatrixXd solutions;
readSolutions(solutions);
cout<<endl<<"Example 4 solutions"<<endl;
cout<<"---------------------"<<endl;
cout << solutions<< endl;
}

TEST_F(MainProgramTests, TestExample5)
{
char example[50] = "../../Examples/FromFileExamples/Example5";
std::cout.setstate(std::ios_base::failbit);
solveProblem(example, 5, 5, 7, false);
std::cout.clear();
MatrixXd solutions;
readSolutions(solutions);
cout<<endl<<"Example 5 solutions"<<endl;
cout<<"---------------------"<<endl;
cout << solutions<< endl;
}
TEST_F(MainProgramTests, TestExample6)
{
char example[50] = "../../Examples/FromFileExamples/Example6";
std::cout.setstate(std::ios_base::failbit);
solveProblem(example, 4, 4, 7, false);
std::cout.clear();
MatrixXd solutions;
readSolutions(solutions);
cout<<endl<<"Example 6 solutions"<<endl;
cout<<"---------------------"<<endl;
cout << solutions<< endl;
}
TEST_F(MainProgramTests, TestExample7)
{
char example[50] = "../../Examples/FromFileExamples/Example7";
std::cout.setstate(std::ios_base::failbit);
solveProblem(example, 5, 4, 7, false);
std::cout.clear();
MatrixXd solutions;
readSolutions(solutions);
cout<<endl<<"Example 7 solutions"<<endl;
cout<<"---------------------"<<endl;
cout << solutions<< endl;
}
TEST_F(MainProgramTests, TestExample8)
{
char example[50] = "../../Examples/FromFileExamples/Example8";
std::cout.setstate(std::ios_base::failbit);
solveProblem(example, 5, 5, 7, false);
std::cout.clear();
MatrixXd solutions;
readSolutions(solutions);
cout<<endl<<"Example 8 solutions"<<endl;
cout<<"---------------------"<<endl;
cout << solutions<< endl;
}


