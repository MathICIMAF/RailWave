#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <complex>
#include <QList>
#include <QFile>
#include <QTextStream>
#include <QDebug>
#include "eigen/linker.h"

using namespace std;

class FComplexMatrix;
class CComplexMatrix;


class FMatrix:public QObject
{
    public:
        int rows, cols;
        double *values;

        FMatrix(QObject*parent=0);
        FMatrix(const FMatrix& mat,QObject*parent=0);
        FMatrix(int prows, int pcols,QObject*parent=0);

        ~FMatrix(){delete[]values;}

        FMatrix transpuesta();
        QList<double> Diagonal();
        FMatrix multi_escalar(double);
        void multi_escalar_this(double);
        void AddThis(const FMatrix &matr2);
        FMatrix operator+ (const FMatrix &matr2);
        FMatrix operator- (const FMatrix &matr2);
        FMatrix operator* (const FMatrix &matr2);
        FMatrix& operator= (const FMatrix& mat);
        void multi_escalar(double esc, FMatrix *result);
        void Add(FComplexMatrix *matr2, FComplexMatrix *result);
        void multi_escalar(complex<double> esc, FComplexMatrix *result);

        void Clear();
        void Symmetrize();
        void AntiSymmetrize();
        complex<double> * ToComplexArray();
        double GetElem(int fila, int columna);
        void ChangeDimensions(int row, int col);
        void SetElem(int fila, int columna, double elem);
        void TrimMatrix(QList<int> zeroBoundaryIndices, FMatrix *result);

        void WriteMatrix(QString path);
};

class CMatrix:public QObject
{
    public:
        double **elems;
        QList<int> neighbors;
        int rows, columns, nrows, ncols;

        CMatrix(int rows, int columns, QList<int> neighbors,QObject*parent=0);

        ~CMatrix() {delete []elems;}

        void multi_escalar(double, CMatrix *result);
        void Add(CComplexMatrix *mtr2, CComplexMatrix *result);
        void multi_escalar(complex<double> val, CComplexMatrix *result);

        void Clear();
        void Symmetrize();
        void AntiSymmetrize();
        QList<int> GetElems(int i0);
        complex<double> * ToComplexArray();
        void SetElems(FMatrix local, int i0, int i1, int i2);
        void TrimMatrix(QList<int> zeroBoundaryIndices, CMatrix *result);
        void SetElems(FMatrix local, int i0, int i1, int i2, int i3, int i4, int i5);

        void WriteMatrix(QString path);
        void WriteMatrixAndNeighbors(QString pathMatrix, QString pathNeigh);
};

class FComplexMatrix:public QObject
{
public:
    int rows, cols;
    complex<double> *values;

    FComplexMatrix(QObject* parent=0);
    FComplexMatrix(int,int,QObject* parent=0);
    FComplexMatrix(const FComplexMatrix& mat,QObject*parent=0);
    FComplexMatrix(complex<double> *values, int rows, int cols,QObject* parent=0);

    ~FComplexMatrix(){delete[] values;}

    void Clear();
    void Hermitize();
    complex<double> * GetColumn(int ith);
    complex<double> GetElem(int fila, int columna);
    void SetElem(int fila, int columna, complex<double> elem);

    void WriteRealComplexMatrix(QString path);
    void WriteImagComplexMatrix(QString path);
};

class CComplexMatrix:public QObject
{
public:
    QList<int> neighbors;
    complex<double> **values;
    int rows, cols, nrows, ncols;

    CComplexMatrix(int rows, int columns, QList<int> neighbors,QObject* parent=0);

    ~CComplexMatrix() {delete values;}

    void Clear();
    complex<double> * ToArray();

    void WriteRealComplexMatrix(QString path);
    void WriteImagComplexMatrix(QString path);
};

int Sign(double x);
QList<int> Sort(QList<int>);

#endif // MATRIX_H
