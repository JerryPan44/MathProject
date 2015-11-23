
class Parser{
    static void readFromFile(char filename[], char*& polynomial1, char*& polynomial2);
    static void readFromStdin(char*& polynomial1, char*& polynomial2);
    static int countLine(FILE * f);
    static void readData(FILE * input, char*& polynomial1, char*& polynomial2);
    static void readLine(FILE *, char *);
public:
    static void readInput(char filename[], char*& polynomial1, char*& polynomial2);

};
