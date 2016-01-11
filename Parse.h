#include "eigen/Eigen/Dense"
#include "eigen/Eigen/Core"
#include "eigen/Eigen/SVD"
#include "eigen/Eigen/Eigenvalues"
class Parser{
    static void readFromFile(char filename[], char*& polynomial1, char*& polynomial2);
    static void readFromStdin(char*& polynomial1, char*& polynomial2);
    static int countLine(FILE * f);
    static void readData(FILE * input, char*& polynomial1, char*& polynomial2);
    static void readLine(FILE *, char *);
public:
    static void readPoints(Eigen::MatrixXd & pointsMatrix, Eigen::MatrixXd & pointsMatrix2);
    static void readInput(char filename[], char*& polynomial1, char*& polynomial2);

};
