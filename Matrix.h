#ifndef _MATRIX_
#define _MATRIX_


class MyMatrix
{
    int Rowdimension;
    int Coldimension;
    double ** matrix;
public:
    MyMatrix(int dim);
    double ** getMatrix()
    {
        return this->matrix;
    }
    MyMatrix(int rowDim, int colDim);
    void Print();
    void generate();
    int getRowDimension()
    {
        return Rowdimension;
    }
    int getColDimension()
    {
        return Coldimension;
    }

};

#endif