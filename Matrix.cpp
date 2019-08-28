#include "Matrix.h"

FMatrix::FMatrix(QObject*parent):QObject(parent)
{
    rows = cols = 0;
    values = new double[0];
}
FMatrix::FMatrix(const FMatrix &mat,QObject*parent):QObject(parent){
    this->cols = mat.cols;
    this->rows = mat.rows;
    long n = rows*cols;
    values = new double[n];
    for(long i = 0; i < n; i++)
        values[i] = mat.values[i];
}

FMatrix::FMatrix(int prows, int pcols, QObject *parent):QObject(parent)
{
    rows = prows;
    cols = pcols;    
    values = new double[rows*cols];
    for(long i = 0 ; i < rows*cols; i ++)
        values[i] = 0;
}

FMatrix& FMatrix::operator =(const FMatrix& mat){
    if(rows != 0 && cols != 0)
        delete[] values;
    this->cols = mat.cols;
    this->rows = mat.rows;
    long n = (mat.rows)*(mat.cols);
    values = new double[n];
    for(long i = 0; i < n; i++)
        this->values[i] = mat.values[i];
    return *this;
}

FMatrix FMatrix::transpuesta()
{
    FMatrix result(cols,rows);
    for(int i = 0; i < cols; i++)
    {
        for(int j = 0; j < rows; j++)
        {
            float elem = GetElem(j,i);
            result.SetElem(i,j,elem);
        }
    }
    return result;
}
QList<double> FMatrix::Diagonal()
{
    QList<double> result;
    for(int i = 0; i < min(rows, cols); i++)
        result.push_back(GetElem(i,i));
    return result;
}
FMatrix FMatrix::multi_escalar(double esc)
{
    FMatrix result(rows,cols);
    for(long i = 0 ; i < rows*cols; i++)
        result.values[i] = values[i]*esc;
    return result;
}
void FMatrix::multi_escalar_this(double esc)
{
    for(long i = 0 ; i < rows*cols;i++)
        values[i] *= esc;
}
void FMatrix::AddThis(const FMatrix &matr2)
{
    for(long i = 0 ; i < rows*cols;i++)
        values[i] += matr2.values[i];
}
FMatrix FMatrix::operator+ (const FMatrix &matr2)
{
    FMatrix result(rows,cols);
    for(long i = 0 ; i < rows*cols;i++)
        result.values[i] = matr2.values[i]+values[i];
    return result;
}
FMatrix FMatrix::operator- (const FMatrix &matr2)
{
    FMatrix result(rows,cols);
    for(int i = 0 ; i < rows*cols;i++)
        result.values[i] = values[i] - matr2.values[i];
    return result;
}
FMatrix FMatrix::operator* (const FMatrix &matr2)
{
    FMatrix result(rows,matr2.cols);
    for(long i = 0; i < rows;i++)
    {
        for(long j = 0; j < matr2.cols; j ++)
        {
            float total = 0;
            for(long k = 0 ; k < cols;k++)
            {
                float val = GetElem(i,k);
                float val2 = matr2.values[k+j*matr2.rows];
                total += val*val2;
            }
            result.SetElem(i,j,total);
        }
    }
    return result;
}
void FMatrix::multi_escalar(double esc, FMatrix *result)
{
    for(long i = 0 ; i < rows*cols;i++)
        result->values[i] = values[i]*esc;
}
void FMatrix::Add(FComplexMatrix *matr2, FComplexMatrix *result)
{
    for(long i = 0 ; i < rows*cols;i++)
    {
        result->values[i].real(matr2->values[i].real() + values[i]);
        result->values[i].imag(matr2->values[i].imag());
    }
}
void FMatrix::multi_escalar(complex<double> esc, FComplexMatrix *result)
{
    for(long i = 0 ; i < rows*cols;i++)
    {
        result->values[i].real(values[i]*esc.real());
        result->values[i].imag(values[i]*esc.imag());
    }
}
void FMatrix::Clear()
{
    for(int i = 0; i < rows*cols; i++)
        values[i] = 0;
}
void FMatrix::Symmetrize()
{
    for(int i = 1; i < rows; i++)
        for(int j = 0; j < i; j++)
            SetElem(i,j,GetElem(j,i));
}
void FMatrix::AntiSymmetrize()
{
    for(int i = 1; i < rows; i++)
        for(int j = 0; j < i; j++)
            SetElem(i,j,-GetElem(j,i));
}
complex<double> * FMatrix::ToComplexArray()
{
    complex<double>* result = new complex<double>[rows*cols];
    for(long i = 0 ; i < rows*cols;i++)
    {
        complex<double> val;
        val.real((double)values[i]);
        result[i] = val;
    }
    return result;
}
double FMatrix::GetElem(int fila, int columna)
{
    return values[fila+columna*rows];
}
void FMatrix::ChangeDimensions(int prows, int pcols)
{
    rows = prows;
    cols = pcols;
    delete(values);
    values = new double[rows*cols];
    for(int i = 0 ; i < rows*cols; i ++)
        values[i] = 0;
}
void FMatrix::SetElem(int fila, int columna, double elem)
{
    values[fila+columna*rows] = elem;
}
void FMatrix::TrimMatrix(QList<int> zeroBoundaryIndices, FMatrix *result)
{
    zeroBoundaryIndices = Sort(zeroBoundaryIndices);    
    int index = 0;
    for(int i = 0 ; i < rows*cols; i++)
    {
        int column = i/rows;
        int row = i%rows;
        if(zeroBoundaryIndices.contains(row) || zeroBoundaryIndices.contains(column))
            result->values[index++] = 0;
        else
            result->values[index++] = values[i];
    }
}
void FMatrix::WriteMatrix(QString path)
{
    FILE *file = fopen(path.toLatin1().constData(),"w");
    for(int i=0;i<rows;i++)
    {
        for(int j=0;j<cols;j++)
        {
            if(GetElem(i,j)>=0)
            {
                if(GetElem(i,j)<10)
                    fprintf(file, "%20.15f ",GetElem(i,j));
                else fprintf(file, "%20.14f ",GetElem(i,j));
            }
            else
            {
                if(GetElem(i,j)>-10)
                    fprintf(file, "%20.14f ",GetElem(i,j));
                else fprintf(file, "%20.13f ",GetElem(i,j));
            }
        }
        fprintf(file,"\n");
    }
    fclose(file);
    //delete file;
}


CMatrix::CMatrix(int prows, int pcolumns, QList<int> pneighbors,QObject*parent):QObject(parent)
{
    rows = prows;
    columns = pcolumns;
    elems = new double*[rows];
    nrows = rows/3;
    ncols = columns/3-1;
    neighbors = pneighbors;

    for(int i = 0 ; i < rows; i ++)
    {
        elems[i] = new double[columns];
        for(int j = 0; j< columns; j++)
            elems[i][j] = 0;
    }
}
void CMatrix::multi_escalar(double val, CMatrix *result)
{
    for(int r = 0; r < rows; r++)
        for(int c = 0; c < columns; c++)
            result->elems[r][c] = elems[r][c]*val;
}
void CMatrix::Add(CComplexMatrix *mtr2, CComplexMatrix *result)
{
    //Check that all three matrices and their neighbors' matrices has the same number of rows and columns
    if(rows != mtr2->rows || rows != result->rows || columns != mtr2->cols || columns != result->cols ||
       nrows != mtr2->nrows || nrows != result->nrows || ncols != mtr2->ncols || ncols != result->ncols)
    {
        qDebug() << "Matrices cannot be added. Dimensions does not agree. \n";
        return;
    }

    //Check that all three matrices has the same structure of neighbors.
    //If not, this operation cannot be done using this method.
    for(int i = 0 ; i < nrows*ncols; i++)
    {
        if(neighbors[i] != mtr2->neighbors[i] || neighbors[i] != result->neighbors[i])
        {
            qDebug()<<"Matrices cannot be added. They has no the same neighbors' structure. \n";
            return;
        }
    }

    for(int r = 0; r < rows; r++)
    {
        for(int c = 0; c < columns; c++)
        {
            result->values[r][c].real(elems[r][c] + mtr2->values[r][c].real());
            result->values[r][c].imag(mtr2->values[r][c].imag());
        }
    }
}
void CMatrix::multi_escalar(complex<double> val, CComplexMatrix *result)
{
    for(int r = 0; r < rows; r++)
    {
        for(int c = 0; c < columns; c++)
        {
            result->values[r][c].real(elems[r][c]*val.real());
            result->values[r][c].imag(elems[r][c]*val.imag());
        }
    }
}
void CMatrix::Clear()
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < columns; j++)
            elems[i][j] = 0;
}
void CMatrix::Symmetrize()
{
    for(int i = 1; i < nrows; i++)
    {
        int threei = 3*i;
        QList<int> neighrs = GetElems(i);
        for(int j = 0; j < neighrs.count(); j++)
        {
            int threej = 3*(j+1);
            if(neighrs[j] == -1)continue;
            //Esta posicion esta por debajo de la diagonal.
            //Es necesario copiar su elemento simetrico.
            if(threei > 3*neighrs[j])
            {
                //obtener la fila a partir de la cual comienzan los
                //datos del "j-esimo vecino de i".
                int threeJvecino = 3*neighrs[j];
                //lista de los vecinos del "j-esimo vecino de i".
                //Lo que se quiere saber es el indice que ocupa i en
                //la lista de los vecinos del "j-esimo vesino de i",
                //para asi poder copiar el pedazo de la matriz
                // (threeJvecino:threeJvecino+2, indIenVecinosJ:indIenVecinosJ+2)
                QList<int> neighrsJneighrs = GetElems(neighrs[j]);
                //Se le suma 1 a indINeighsJ proque debemos recordar
                //que en las columnas 0:2 de la matriz compactas se
                //almacena la relacion del vertice consigo mismo.
                int threeIndINeighsJ = 3*(neighrsJneighrs.indexOf(i) + 1);
                elems[threei][threej] = elems[threeJvecino][threeIndINeighsJ];
                elems[threei+1][threej] = elems[threeJvecino][threeIndINeighsJ+1];
                elems[threei+2][threej] = elems[threeJvecino][threeIndINeighsJ+2];

                elems[threei][threej+1] = elems[threeJvecino+1][threeIndINeighsJ];
                elems[threei+1][threej+1] = elems[threeJvecino+1][threeIndINeighsJ+1];
                elems[threei+2][threej+1] = elems[threeJvecino+1][threeIndINeighsJ+2];

                elems[threei][threej+2] = elems[threeJvecino+2][threeIndINeighsJ];
                elems[threei+1][threej+2] = elems[threeJvecino+2][threeIndINeighsJ+1];
                elems[threei+2][threej+2] = elems[threeJvecino+2][threeIndINeighsJ+2];
            }
        }
    }
}
void CMatrix::AntiSymmetrize()
{
    for(int i = 1; i < nrows; i++)
    {
        int threei = 3*i;
        QList<int> neighrs = GetElems(i);
        for(int j = 0; j < neighrs.count(); j++)
        {
            int threej = 3*(j+1);
            if(neighrs[j] == -1)continue;
            //Esta posicion esta por debajo de la diagonal.
            //Es necesario copiar su elemento simetrico.
            if(threei > 3*neighrs[j])
            {
                //obtener la fila a partir de la cual comienzan los
                //datos del "j-esimo vecino de i".
                int threeJvecino = 3*neighrs[j];
                //lista de los vecinos del "j-esimo vecino de i".
                //Lo que se quiere saber es el indice que ocupa i en
                //la lista de los vecinos del "j-esimo vesino de i",
                //para asi poder copiar el pedazo de la matriz
                // (threeJvecino:threeJvecino+2, indIenVecinosJ:indIenVecinosJ+2)
                QList<int> neighrsJneighrs = GetElems(neighrs[j]);
                //Se le suma 1 a indINeighsJ proque debemos recordar
                //que en las columnas 0:2 de la matriz compactas se
                //almacena la relacion del vertice consigo mismo.
                int threeIndINeighsJ = 3*(neighrsJneighrs.indexOf(i) + 1);

                elems[threei][threej] = -elems[threeJvecino][threeIndINeighsJ];
                elems[threei+1][threej] = -elems[threeJvecino][threeIndINeighsJ+1];
                elems[threei+2][threej] = -elems[threeJvecino][threeIndINeighsJ+2];

                elems[threei][threej+1] = -elems[threeJvecino+1][threeIndINeighsJ];
                elems[threei+1][threej+1] = -elems[threeJvecino+1][threeIndINeighsJ+1];
                elems[threei+2][threej+1] = -elems[threeJvecino+1][threeIndINeighsJ+2];

                elems[threei][threej+2] = -elems[threeJvecino+2][threeIndINeighsJ];
                elems[threei+1][threej+2] = -elems[threeJvecino+2][threeIndINeighsJ+1];
                elems[threei+2][threej+2] = -elems[threeJvecino+2][threeIndINeighsJ+2];
            }
        }
    }
}
QList<int> CMatrix::GetElems(int i0)
{
    QList<int> toreturn;
    for(int i = 0 ; i < ncols; i++)
        toreturn.push_back(neighbors[nrows*i+i0]);
    return toreturn;
}
complex<double> * CMatrix::ToComplexArray()
{
    complex<double> *result = new complex<double>[rows*columns];
    int index = 0;
    for(int i = 0 ; i < columns; i++)
    {
        for(int j = 0; j < rows; j ++)
        {
            result[index].real(elems[j][i]);
            index++;
        }
    }
    return result;
}
void CMatrix::SetElems(FMatrix local, int i0, int i1, int i2)
{
    /*
     Este metodo llena los valores en la matriz compacta
     correspondientes a los nodos i0, i1 e i2. Se van a
     guardar solo las posiciones que se encuentren por encima
     de la diagonal. Para obtener la matriz simetrica hay que llamar
     al método Symmetrize.
    */
    int threei0 = 3*i0;
    elems[threei0][0] = elems[threei0][0]+ local.GetElem(0,0);
    elems[threei0+1][0] = elems[threei0+1][0]+ local.GetElem(1,0);
    elems[threei0+2][0] = elems[threei0+2][0]+ local.GetElem(2,0);
    elems[threei0][1] = elems[threei0][1]+ local.GetElem(0,1);
    elems[threei0+1][1] = elems[threei0+1][1]+ local.GetElem(1,1);
    elems[threei0+2][1] = elems[threei0+2][1]+ local.GetElem(2,1);
    elems[threei0][2] = elems[threei0][2]+ local.GetElem(0,2);
    elems[threei0+1][2] = elems[threei0+1][2]+ local.GetElem(1,2);
    elems[threei0+2][2] = elems[threei0+2][2]+ local.GetElem(2,2);

    QList<int> neighrs = GetElems(i0);
    int ind1_i0 = neighrs.indexOf(i1);
    int threeind1_i0 = 3*(ind1_i0+1);
    if(threei0 < 3*i1)
    {
        elems[threei0][threeind1_i0] = elems[threei0][threeind1_i0]+ local.GetElem(0,3);
        elems[threei0+1][threeind1_i0] = elems[threei0+1][threeind1_i0]+ local.GetElem(1,3);
        elems[threei0+2][threeind1_i0] = elems[threei0+2][threeind1_i0]+ local.GetElem(2,3);
        elems[threei0][threeind1_i0+1] = elems[threei0][threeind1_i0+1]+ local.GetElem(0,4);
        elems[threei0+1][threeind1_i0+1] = elems[threei0+1][threeind1_i0+1]+ local.GetElem(1,4);
        elems[threei0+2][threeind1_i0+1] = elems[threei0+2][threeind1_i0+1]+ local.GetElem(2,4);
        elems[threei0][threeind1_i0+2] = elems[threei0][threeind1_i0+2]+ local.GetElem(0,5);
        elems[threei0+1][threeind1_i0+2] = elems[threei0+1][threeind1_i0+2]+ local.GetElem(1,5);
        elems[threei0+2][threeind1_i0+2] = elems[threei0+2][threeind1_i0+2]+ local.GetElem(2,5);
    }

    int ind2_i0 = neighrs.indexOf(i2);
    int threeind2_i0 = 3*(ind2_i0+1);
    if(threei0 < 3*i2)
    {
        elems[threei0][threeind2_i0] = elems[threei0][threeind2_i0]+ local.GetElem(0,6);
        elems[threei0+1][threeind2_i0] = elems[threei0+1][threeind2_i0]+ local.GetElem(1,6);
        elems[threei0+2][threeind2_i0] = elems[threei0+2][threeind2_i0]+ local.GetElem(2,6);
        elems[threei0][threeind2_i0+1] = elems[threei0][threeind2_i0+1]+ local.GetElem(0,7);
        elems[threei0+1][threeind2_i0+1] = elems[threei0+1][threeind2_i0+1]+ local.GetElem(1,7);
        elems[threei0+2][threeind2_i0+1] = elems[threei0+2][threeind2_i0+1]+ local.GetElem(2,7);
        elems[threei0][threeind2_i0+2] = elems[threei0][threeind2_i0+2]+ local.GetElem(0,8);
        elems[threei0+1][threeind2_i0+2] = elems[threei0+1][threeind2_i0+2]+ local.GetElem(1,8);
        elems[threei0+2][threeind2_i0+2] = elems[threei0+2][threeind2_i0+2]+ local.GetElem(2,8);
    }


    int threei1 = 3*i1;
    elems[threei1][0] = elems[threei1][0]+ local.GetElem(3,3);
    elems[threei1+1][0] = elems[threei1+1][0]+ local.GetElem(4,3);
    elems[threei1+2][0] = elems[threei1+2][0]+ local.GetElem(5,3);
    elems[threei1][1] = elems[threei1][1]+ local.GetElem(3,4);
    elems[threei1+1][1] = elems[threei1+1][1]+ local.GetElem(4,4);
    elems[threei1+2][1] = elems[threei1+2][1]+ local.GetElem(5,4);
    elems[threei1][2] = elems[threei1][2]+ local.GetElem(3,5);
    elems[threei1+1][2] = elems[threei1+1][2]+ local.GetElem(4,5);
    elems[threei1+2][2] = elems[threei1+2][2]+ local.GetElem(5,5);

    neighrs = GetElems(i1);
    int ind0_i1 = neighrs.indexOf(i0);
    int threeind0_i1 = 3*(ind0_i1+1);
    if(threei1 < 3*i0)
    {
        elems[threei1][threeind0_i1] = elems[threei1][threeind0_i1]+ local.GetElem(3,0);
        elems[threei1+1][threeind0_i1] = elems[threei1+1][threeind0_i1]+ local.GetElem(4,0);
        elems[threei1+2][threeind0_i1] = elems[threei1+2][threeind0_i1]+ local.GetElem(5,0);
        elems[threei1][threeind0_i1+1] = elems[threei1][threeind0_i1+1]+ local.GetElem(3,1);
        elems[threei1+1][threeind0_i1+1] = elems[threei1+1][threeind0_i1+1]+ local.GetElem(4,1);
        elems[threei1+2][threeind0_i1+1] = elems[threei1+2][threeind0_i1+1]+ local.GetElem(5,1);
        elems[threei1][threeind0_i1+2] = elems[threei1][threeind0_i1+2]+ local.GetElem(3,2);
        elems[threei1+1][threeind0_i1+2] = elems[threei1+1][threeind0_i1+2]+ local.GetElem(4,2);
        elems[threei1+2][threeind0_i1+2] = elems[threei1+2][threeind0_i1+2]+ local.GetElem(5,2);
    }

    int ind2_i1 = neighrs.indexOf(i2);
    int threeind2_i1 = 3*(ind2_i1+1);
    if(threei1 < 3*i2)
    {
        elems[threei1][threeind2_i1] = elems[threei1][threeind2_i1]+ local.GetElem(3,6);
        elems[threei1+1][threeind2_i1] = elems[threei1+1][threeind2_i1]+ local.GetElem(4,6);
        elems[threei1+2][threeind2_i1] = elems[threei1+2][threeind2_i1]+ local.GetElem(5,6);
        elems[threei1][threeind2_i1+1] = elems[threei1][threeind2_i1+1]+ local.GetElem(3,7);
        elems[threei1+1][threeind2_i1+1] = elems[threei1+1][threeind2_i1+1]+ local.GetElem(4,7);
        elems[threei1+2][threeind2_i1+1] = elems[threei1+2][threeind2_i1+1]+ local.GetElem(5,7);
        elems[threei1][threeind2_i1+2] = elems[threei1][threeind2_i1+2]+ local.GetElem(3,8);
        elems[threei1+1][threeind2_i1+2] = elems[threei1+1][threeind2_i1+2]+ local.GetElem(4,8);
        elems[threei1+2][threeind2_i1+2] = elems[threei1+2][threeind2_i1+2]+ local.GetElem(5,8);
    }


    int threei2 = 3*i2;
    elems[threei2][0] = elems[threei2][0]+ local.GetElem(6,6);
    elems[threei2+1][0] = elems[threei2+1][0]+ local.GetElem(7,6);
    elems[threei2+2][0] = elems[threei2+2][0]+ local.GetElem(8,6);
    elems[threei2][1] = elems[threei2][1]+ local.GetElem(6,7);
    elems[threei2+1][1] = elems[threei2+1][1]+ local.GetElem(7,7);
    elems[threei2+2][1] = elems[threei2+2][1]+ local.GetElem(8,7);
    elems[threei2][2] = elems[threei2][2]+ local.GetElem(6,8);
    elems[threei2+1][2] = elems[threei2+1][2]+ local.GetElem(7,8);
    elems[threei2+2][2] = elems[threei2+2][2]+ local.GetElem(8,8);

    neighrs = GetElems(i2);
    int ind0_i2 = neighrs.indexOf(i0);
    int threeind0_i2 = 3*(ind0_i2+1);
    if(threei2 < 3*i0)
    {
        elems[threei2][threeind0_i2] = elems[threei2][threeind0_i2]+ local.GetElem(6,0);
        elems[threei2+1][threeind0_i2] = elems[threei2+1][threeind0_i2]+ local.GetElem(7,0);
        elems[threei2+2][threeind0_i2] = elems[threei2+2][threeind0_i2]+ local.GetElem(8,0);
        elems[threei2][threeind0_i2+1] = elems[threei2][threeind0_i2+1]+ local.GetElem(6,1);
        elems[threei2+1][threeind0_i2+1] = elems[threei2+1][threeind0_i2+1]+ local.GetElem(7,1);
        elems[threei2+2][threeind0_i2+1] = elems[threei2+2][threeind0_i2+1]+ local.GetElem(8,1);
        elems[threei2][threeind0_i2+2] = elems[threei2][threeind0_i2+2]+ local.GetElem(6,2);
        elems[threei2+1][threeind0_i2+2] = elems[threei2+1][threeind0_i2+2]+ local.GetElem(7,2);
        elems[threei2+2][threeind0_i2+2] = elems[threei2+2][threeind0_i2+2]+ local.GetElem(8,2);
    }

    int ind1_i2 = neighrs.indexOf(i1);
    int threeind1_i2 = 3*(ind1_i2+1);
    if(threei2 < 3*i1)
    {
        elems[threei2][threeind1_i2] = elems[threei2][threeind1_i2]+ local.GetElem(6,3);
        elems[threei2+1][threeind1_i2] = elems[threei2+1][threeind1_i2]+ local.GetElem(7,3);
        elems[threei2+2][threeind1_i2] = elems[threei2+2][threeind1_i2]+ local.GetElem(8,3);
        elems[threei2][threeind1_i2+1] = elems[threei2][threeind1_i2+1]+ local.GetElem(6,4);
        elems[threei2+1][threeind1_i2+1] = elems[threei2+1][threeind1_i2+1]+ local.GetElem(7,4);
        elems[threei2+2][threeind1_i2+1] = elems[threei2+2][threeind1_i2+1]+ local.GetElem(8,4);
        elems[threei2][threeind1_i2+2] = elems[threei2][threeind1_i2+2]+ local.GetElem(6,5);
        elems[threei2+1][threeind1_i2+2] = elems[threei2+1][threeind1_i2+2]+ local.GetElem(7,5);
        elems[threei2+2][threeind1_i2+2] = elems[threei2+2][threeind1_i2+2]+ local.GetElem(8,5);
    }
}
void CMatrix::TrimMatrix(QList<int> zeroBoundaryIndices, CMatrix *result)
{
    zeroBoundaryIndices = Sort(zeroBoundaryIndices);
    for(int i = 0; i < nrows; i++)
    {
        if(zeroBoundaryIndices.contains(3*i))
        {
            for(int j = 0; j < columns; j++)
                result->elems[3*i][j] = 0;
        }
        else
        {
            result->elems[3*i][3*i] = elems[3*i][3*i];
            result->elems[3*i][3*i+1] = elems[3*i][3*i+1];
            result->elems[3*i][3*i+2] = elems[3*i][3*i+2];

            result->elems[3*i+1][3*i] = elems[3*i+1][3*i];
            result->elems[3*i+1][3*i+1] = elems[3*i+1][3*i+1];
            result->elems[3*i+1][3*i+2] = elems[3*i+1][3*i+2];

            result->elems[3*i+2][3*i] = elems[3*i+2][3*i];
            result->elems[3*i+2][3*i+1] = elems[3*i+2][3*i+1];
            result->elems[3*i+2][3*i+2] = elems[3*i+2][3*i+2];

            for(int j = 0; j < ncols  ; j++)
            {
                int indexNeighbor = neighbors[(j*nrows) + i];
                if(indexNeighbor == -1)break;
                int c = 3*indexNeighbor;
                //AQUI se le suma 1 a j (3*j+1) porque en la matriz compacta, el j-esimo vecino esta en las columnas 3*(j+1),
                //3*(j+1)+1 y 3*(j+1)+2 ya que en las 3 primeras columnas esta la relacion del vertices con el mismo.
                if(zeroBoundaryIndices.contains(c))
                    result->elems[3*i][3*(j+1)] = 0;
                else
                    result->elems[3*i][3*(j+1)] = elems[3*i][3*(j+1)];

                if(zeroBoundaryIndices.contains(c+1))
                    result->elems[3*i][3*(j+1)+1] = 0;
                else
                    result->elems[3*i][3*(j+1)+1] = elems[3*i][3*(j+1)+1];

                 if(zeroBoundaryIndices.contains(c+2))
                     result->elems[3*i][3*(j+1)+2] = 0;
                 else result->elems[3*i][3*(j+1)+2] = elems[3*i][3*(j+1)+2];
             }
        }
        if(zeroBoundaryIndices.contains(3*i+1))
        {
            for(int j = 0; j < columns; j++)
                result->elems[3*i+1][j] = 0;
        }
        else
        {
            for(int j = 0; j < ncols; j++)
            {
                int indexNeighbor = neighbors[(j*i) + i];
                int c = 3*indexNeighbor;
                //AQUI se le suma 1 a j (3*j+1) porque en la matriz compacta, el j-esimo vecino esta en las columnas 3*(j+1),
                //3*(j+1)+1 y 3*(j+1)+2 ya que en las 3 primeras columnas esta la relacion del vertices con el mismo.
                if(zeroBoundaryIndices.contains(c))
                    result->elems[3*i+1][3*(j+1)] = 0;
                else
                    result->elems[3*i+1][3*(j+1)] = elems[3*i+1][3*(j+1)];

                if(zeroBoundaryIndices.contains(c+1))
                    result->elems[3*i+1][3*(j+1)+1] = 0;
                else
                    result->elems[3*i+1][3*(j+1)+1] = elems[3*i+1][3*(j+1)+1];

                 if(zeroBoundaryIndices.contains(c+2))
                     result->elems[3*i+1][3*(j+1)+2] = 0;
                 else result->elems[3*i+1][3*(j+1)+2] = elems[3*i+1][3*(j+1)+2];
             }
        }
        if(zeroBoundaryIndices.contains(3*i+2))
        {
            for(int j = 0; j < columns; j++)
                result->elems[3*i+2][j] = 0;
        }
        else
        {
            for(int j = 0; j < ncols; j++)
            {
                int indexNeighbor = neighbors[(j*i) + i];
                int c = 3*indexNeighbor;
                //AQUI se le suma 1 a j (3*j+1) porque en la matriz compacta, el j-esimo vecino esta en las columnas 3*(j+1),
                //3*(j+1)+1 y 3*(j+1)+2 ya que en las 3 primeras columnas esta la relacion del vertices con el mismo.
                if(zeroBoundaryIndices.contains(c))
                    result->elems[3*i+2][3*(j+1)] = 0;
                else
                    result->elems[3*i+2][3*(j+1)] = elems[3*i+2][3*(j+1)];

                if(zeroBoundaryIndices.contains(c+1))
                    result->elems[3*i+2][3*(j+1)+1] = 0;
                else
                    result->elems[3*i+2][3*(j+1)+1] = elems[3*i+2][3*(j+1)+1];

                 if(zeroBoundaryIndices.contains(c+2))
                     result->elems[3*i+2][3*(j+1)+2] = 0;
                 else result->elems[3*i+2][3*(j+1)+2] = elems[3*i+2][3*(j+1)+2];
             }
        }
    }
}
void CMatrix::SetElems(FMatrix local, int i0, int i1, int i2, int i3, int i4, int i5)
{
    /*
         Este metodo llena los valores en la matriz compacta
         correspondientes a los nodos i0, i1, i2, i3, i4 e i5. Se van a
         guardar solo las posiciones que se encuentren por encima
         de la diagonal. Para obtener la matriz simetrica hay que llamar
         al método Symmetrize.
        */
    ///////////// Llenar las posiciones de la matriz correspondientes al vertice i0 /////////////
    int threei0 = 3*i0;
    elems[threei0][0] = elems[threei0][0]+ local.GetElem(0,0);
    elems[threei0+1][0] = elems[threei0+1][0]+ local.GetElem(1,0);
    elems[threei0+2][0] = elems[threei0+2][0]+ local.GetElem(2,0);
    elems[threei0][1] = elems[threei0][1]+ local.GetElem(0,1);
    elems[threei0+1][1] = elems[threei0+1][1]+ local.GetElem(1,1);
    elems[threei0+2][1] = elems[threei0+2][1]+ local.GetElem(2,1);
    elems[threei0][2] = elems[threei0][2]+ local.GetElem(0,2);
    elems[threei0+1][2] = elems[threei0+1][2]+ local.GetElem(1,2);
    elems[threei0+2][2] = elems[threei0+2][2]+ local.GetElem(2,2);

    QList<int> neighrs = GetElems(i0);
    int ind1_i0 = neighrs.indexOf(i1);
    if(threei0 < 3*i1)
    {
        int threeind1_i0 = 3*(ind1_i0+1);
        elems[threei0][threeind1_i0] = elems[threei0][threeind1_i0]+ local.GetElem(0,3);
        elems[threei0+1][threeind1_i0] = elems[threei0+1][threeind1_i0]+ local.GetElem(1,3);
        elems[threei0+2][threeind1_i0] = elems[threei0+2][threeind1_i0]+ local.GetElem(2,3);
        elems[threei0][threeind1_i0+1] = elems[threei0][threeind1_i0+1]+ local.GetElem(0,4);
        elems[threei0+1][threeind1_i0+1] = elems[threei0+1][threeind1_i0+1]+ local.GetElem(1,4);
        elems[threei0+2][threeind1_i0+1] = elems[threei0+2][threeind1_i0+1]+ local.GetElem(2,4);
        elems[threei0][threeind1_i0+2] = elems[threei0][threeind1_i0+2]+ local.GetElem(0,5);
        elems[threei0+1][threeind1_i0+2] = elems[threei0+1][threeind1_i0+2]+ local.GetElem(1,5);
        elems[threei0+2][threeind1_i0+2] = elems[threei0+2][threeind1_i0+2]+ local.GetElem(2,5);
    }

    int ind2_i0 = neighrs.indexOf(i2);
    if(threei0 < 3*i2)
    {
        int threeind2_i0 = 3*(ind2_i0+1);
        elems[threei0][threeind2_i0] = elems[threei0][threeind2_i0]+ local.GetElem(0,6);
        elems[threei0+1][threeind2_i0] = elems[threei0+1][threeind2_i0]+ local.GetElem(1,6);
        elems[threei0+2][threeind2_i0] = elems[threei0+2][threeind2_i0]+ local.GetElem(2,6);
        elems[threei0][threeind2_i0+1] = elems[threei0][threeind2_i0+1]+ local.GetElem(0,7);
        elems[threei0+1][threeind2_i0+1] = elems[threei0+1][threeind2_i0+1]+ local.GetElem(1,7);
        elems[threei0+2][threeind2_i0+1] = elems[threei0+2][threeind2_i0+1]+ local.GetElem(2,7);
        elems[threei0][threeind2_i0+2] = elems[threei0][threeind2_i0+2]+ local.GetElem(0,8);
        elems[threei0+1][threeind2_i0+2] = elems[threei0+1][threeind2_i0+2]+ local.GetElem(1,8);
        elems[threei0+2][threeind2_i0+2] = elems[threei0+2][threeind2_i0+2]+ local.GetElem(2,8);
    }

    int ind3_i0 = neighrs.indexOf(i3);
    if(threei0 < 3*i3)
    {
        int threeind3_i0 = 3*(ind3_i0+1);
        elems[threei0][threeind3_i0] = elems[threei0][threeind3_i0]+ local.GetElem(0,9);
        elems[threei0+1][threeind3_i0] = elems[threei0+1][threeind3_i0]+ local.GetElem(1,9);
        elems[threei0+2][threeind3_i0] = elems[threei0+2][threeind3_i0]+ local.GetElem(2,9);
        elems[threei0][threeind3_i0+1] = elems[threei0][threeind3_i0+1]+ local.GetElem(0,10);
        elems[threei0+1][threeind3_i0+1] = elems[threei0+1][threeind3_i0+1]+ local.GetElem(1,10);
        elems[threei0+2][threeind3_i0+1] = elems[threei0+2][threeind3_i0+1]+ local.GetElem(2,10);
        elems[threei0][threeind3_i0+2] = elems[threei0][threeind3_i0+2]+ local.GetElem(0,11);
        elems[threei0+1][threeind3_i0+2] = elems[threei0+1][threeind3_i0+2]+ local.GetElem(1,11);
        elems[threei0+2][threeind3_i0+2] = elems[threei0+2][threeind3_i0+2]+ local.GetElem(2,11);
    }

    int ind4_i0 = neighrs.indexOf(i4);
    if(threei0 < 3*i4)
    {
        int threeind4_i0 = 3*(ind4_i0+1);
        elems[threei0][threeind4_i0] = elems[threei0][threeind4_i0]+ local.GetElem(0,12);
        elems[threei0+1][threeind4_i0] = elems[threei0+1][threeind4_i0]+ local.GetElem(1,12);
        elems[threei0+2][threeind4_i0] = elems[threei0+2][threeind4_i0]+ local.GetElem(2,12);
        elems[threei0][threeind4_i0+1] = elems[threei0][threeind4_i0+1]+ local.GetElem(0,13);
        elems[threei0+1][threeind4_i0+1] = elems[threei0+1][threeind4_i0+1]+ local.GetElem(1,13);
        elems[threei0+2][threeind4_i0+1] = elems[threei0+2][threeind4_i0+1]+ local.GetElem(2,13);
        elems[threei0][threeind4_i0+2] = elems[threei0][threeind4_i0+2]+ local.GetElem(0,14);
        elems[threei0+1][threeind4_i0+2] = elems[threei0+1][threeind4_i0+2]+ local.GetElem(1,14);
        elems[threei0+2][threeind4_i0+2] = elems[threei0+2][threeind4_i0+2]+ local.GetElem(2,14);
    }

    int ind5_i0 = neighrs.indexOf(i5);
    if(threei0 < 3*i5)
    {
        int threeind5_i0 = 3*(ind5_i0+1);
        elems[threei0][threeind5_i0] = elems[threei0][threeind5_i0]+ local.GetElem(0,15);
        elems[threei0+1][threeind5_i0] = elems[threei0+1][threeind5_i0]+ local.GetElem(1,15);
        elems[threei0+2][threeind5_i0] = elems[threei0+2][threeind5_i0]+ local.GetElem(2,15);
        elems[threei0][threeind5_i0+1] = elems[threei0][threeind5_i0+1]+ local.GetElem(0,16);
        elems[threei0+1][threeind5_i0+1] = elems[threei0+1][threeind5_i0+1]+ local.GetElem(1,16);
        elems[threei0+2][threeind5_i0+1] = elems[threei0+2][threeind5_i0+1]+ local.GetElem(2,16);
        elems[threei0][threeind5_i0+2] = elems[threei0][threeind5_i0+2]+ local.GetElem(0,17);
        elems[threei0+1][threeind5_i0+2] = elems[threei0+1][threeind5_i0+2]+ local.GetElem(1,17);
        elems[threei0+2][threeind5_i0+2] = elems[threei0+2][threeind5_i0+2]+ local.GetElem(2,17);
    }

    ///////////// Llenar las posiciones de la matriz correspondientes al vertice i1 /////////////
    int threei1 = 3*i1;
    elems[threei1][0] = elems[threei1][0]+ local.GetElem(3,3);
    elems[threei1+1][0] = elems[threei1+1][0]+ local.GetElem(4,3);
    elems[threei1+2][0] = elems[threei1+2][0]+ local.GetElem(5,3);
    elems[threei1][1] = elems[threei1][1]+ local.GetElem(3,4);
    elems[threei1+1][1] = elems[threei1+1][1]+ local.GetElem(4,4);
    elems[threei1+2][1] = elems[threei1+2][1]+ local.GetElem(5,4);
    elems[threei1][2] = elems[threei1][2]+ local.GetElem(3,5);
    elems[threei1+1][2] = elems[threei1+1][2]+ local.GetElem(4,5);
    elems[threei1+2][2] = elems[threei1+2][2]+ local.GetElem(5,5);

    neighrs = GetElems(i1);
    int ind0_i1 = neighrs.indexOf(i0);
    if(threei1 < 3*i0)
    {
        int threeind0_i1 = 3*(ind0_i1+1);
        elems[threei1][threeind0_i1] = elems[threei1][threeind0_i1]+ local.GetElem(3,0);
        elems[threei1+1][threeind0_i1] = elems[threei1+1][threeind0_i1]+ local.GetElem(4,0);
        elems[threei1+2][threeind0_i1] = elems[threei1+2][threeind0_i1]+ local.GetElem(5,0);
        elems[threei1][threeind0_i1+1] = elems[threei1][threeind0_i1+1]+ local.GetElem(3,1);
        elems[threei1+1][threeind0_i1+1] = elems[threei1+1][threeind0_i1+1]+ local.GetElem(4,1);
        elems[threei1+2][threeind0_i1+1] = elems[threei1+2][threeind0_i1+1]+ local.GetElem(5,1);
        elems[threei1][threeind0_i1+2] = elems[threei1][threeind0_i1+2]+ local.GetElem(3,2);
        elems[threei1+1][threeind0_i1+2] = elems[threei1+1][threeind0_i1+2]+ local.GetElem(4,2);
        elems[threei1+2][threeind0_i1+2] = elems[threei1+2][threeind0_i1+2]+ local.GetElem(5,2);
    }

    int ind2_i1 = neighrs.indexOf(i2);
    if(threei1 < 3*i2)
    {
        int threeind2_i1 = 3*(ind2_i1+1);
        elems[threei1][threeind2_i1] = elems[threei1][threeind2_i1]+ local.GetElem(3,6);
        elems[threei1+1][threeind2_i1] = elems[threei1+1][threeind2_i1]+ local.GetElem(4,6);
        elems[threei1+2][threeind2_i1] = elems[threei1+2][threeind2_i1]+ local.GetElem(5,6);
        elems[threei1][threeind2_i1+1] = elems[threei1][threeind2_i1+1]+ local.GetElem(3,7);
        elems[threei1+1][threeind2_i1+1] = elems[threei1+1][threeind2_i1+1]+ local.GetElem(4,7);
        elems[threei1+2][threeind2_i1+1] = elems[threei1+2][threeind2_i1+1]+ local.GetElem(5,7);
        elems[threei1][threeind2_i1+2] = elems[threei1][threeind2_i1+2]+ local.GetElem(3,8);
        elems[threei1+1][threeind2_i1+2] = elems[threei1+1][threeind2_i1+2]+ local.GetElem(4,8);
        elems[threei1+2][threeind2_i1+2] = elems[threei1+2][threeind2_i1+2]+ local.GetElem(5,8);
    }

    int ind3_i1 = neighrs.indexOf(i3);
    if(threei1 < 3*i3)
    {
        int threeind3_i1 = 3*(ind3_i1+1);
        elems[threei1][threeind3_i1] = elems[threei1][threeind3_i1]+ local.GetElem(3,9);
        elems[threei1+1][threeind3_i1] = elems[threei1+1][threeind3_i1]+ local.GetElem(4,9);
        elems[threei1+2][threeind3_i1] = elems[threei1+2][threeind3_i1]+ local.GetElem(5,9);
        elems[threei1][threeind3_i1+1] = elems[threei1][threeind3_i1+1]+ local.GetElem(3,10);
        elems[threei1+1][threeind3_i1+1] = elems[threei1+1][threeind3_i1+1]+ local.GetElem(4,10);
        elems[threei1+2][threeind3_i1+1] = elems[threei1+2][threeind3_i1+1]+ local.GetElem(5,10);
        elems[threei1][threeind3_i1+2] = elems[threei1][threeind3_i1+2]+ local.GetElem(3,11);
        elems[threei1+1][threeind3_i1+2] = elems[threei1+1][threeind3_i1+2]+ local.GetElem(4,11);
        elems[threei1+2][threeind3_i1+2] = elems[threei1+2][threeind3_i1+2]+ local.GetElem(5,11);
    }

    int ind4_i1 = neighrs.indexOf(i4);
    if(threei1 < 3*i4)
    {
        int threeind4_i1 = 3*(ind4_i1+1);
        elems[threei1][threeind4_i1] = elems[threei1][threeind4_i1]+ local.GetElem(3,12);
        elems[threei1+1][threeind4_i1] = elems[threei1+1][threeind4_i1]+ local.GetElem(4,12);
        elems[threei1+2][threeind4_i1] = elems[threei1+2][threeind4_i1]+ local.GetElem(5,12);
        elems[threei1][threeind4_i1+1] = elems[threei1][threeind4_i1+1]+ local.GetElem(3,13);
        elems[threei1+1][threeind4_i1+1] = elems[threei1+1][threeind4_i1+1]+ local.GetElem(4,13);
        elems[threei1+2][threeind4_i1+1] = elems[threei1+2][threeind4_i1+1]+ local.GetElem(5,13);
        elems[threei1][threeind4_i1+2] = elems[threei1][threeind4_i1+2]+ local.GetElem(3,14);
        elems[threei1+1][threeind4_i1+2] = elems[threei1+1][threeind4_i1+2]+ local.GetElem(4,14);
        elems[threei1+2][threeind4_i1+2] = elems[threei1+2][threeind4_i1+2]+ local.GetElem(5,14);
    }

    int ind5_i1 = neighrs.indexOf(i5);
    if(threei1 < 3*i5)
    {
        int threeind5_i1 = 3*(ind5_i1+1);
        elems[threei1][threeind5_i1] = elems[threei1][threeind5_i1]+ local.GetElem(3,15);
        elems[threei1+1][threeind5_i1] = elems[threei1+1][threeind5_i1]+ local.GetElem(4,15);
        elems[threei1+2][threeind5_i1] = elems[threei1+2][threeind5_i1]+ local.GetElem(5,15);
        elems[threei1][threeind5_i1+1] = elems[threei1][threeind5_i1+1]+ local.GetElem(3,16);
        elems[threei1+1][threeind5_i1+1] = elems[threei1+1][threeind5_i1+1]+ local.GetElem(4,16);
        elems[threei1+2][threeind5_i1+1] = elems[threei1+2][threeind5_i1+1]+ local.GetElem(5,16);
        elems[threei1][threeind5_i1+2] = elems[threei1][threeind5_i1+2]+ local.GetElem(3,17);
        elems[threei1+1][threeind5_i1+2] = elems[threei1+1][threeind5_i1+2]+ local.GetElem(4,17);
        elems[threei1+2][threeind5_i1+2] = elems[threei1+2][threeind5_i1+2]+ local.GetElem(5,17);
    }

    ///////////// Llenar las posiciones de la matriz correspondientes al vertice i2 /////////////
    int threei2 = 3*i2;
    elems[threei2][0] = elems[threei2][0]+ local.GetElem(6,6);
    elems[threei2+1][0] = elems[threei2+1][0]+ local.GetElem(7,6);
    elems[threei2+2][0] = elems[threei2+2][0]+ local.GetElem(8,6);
    elems[threei2][1] = elems[threei2][1]+ local.GetElem(6,7);
    elems[threei2+1][1] = elems[threei2+1][1]+ local.GetElem(7,7);
    elems[threei2+2][1] = elems[threei2+2][1]+ local.GetElem(8,7);
    elems[threei2][2] = elems[threei2][2]+ local.GetElem(6,8);
    elems[threei2+1][2] = elems[threei2+1][2]+ local.GetElem(7,8);
    elems[threei2+2][2] = elems[threei2+2][2]+ local.GetElem(8,8);

    neighrs = GetElems(i2);
    int ind0_i2 = neighrs.indexOf(i0);
    if(threei2 < 3*i0)
    {
        int threeind0_i2 = 3*(ind0_i2+1);
        elems[threei2][threeind0_i2] = elems[threei2][threeind0_i2]+ local.GetElem(6,0);
        elems[threei2+1][threeind0_i2] = elems[threei2+1][threeind0_i2]+ local.GetElem(7,0);
        elems[threei2+2][threeind0_i2] = elems[threei2+2][threeind0_i2]+ local.GetElem(8,0);
        elems[threei2][threeind0_i2+1] = elems[threei2][threeind0_i2+1]+ local.GetElem(6,1);
        elems[threei2+1][threeind0_i2+1] = elems[threei2+1][threeind0_i2+1]+ local.GetElem(7,1);
        elems[threei2+2][threeind0_i2+1] = elems[threei2+2][threeind0_i2+1]+ local.GetElem(8,1);
        elems[threei2][threeind0_i2+2] = elems[threei2][threeind0_i2+2]+ local.GetElem(6,2);
        elems[threei2+1][threeind0_i2+2] = elems[threei2+1][threeind0_i2+2]+ local.GetElem(7,2);
        elems[threei2+2][threeind0_i2+2] = elems[threei2+2][threeind0_i2+2]+ local.GetElem(8,2);
    }

    int ind1_i2 = neighrs.indexOf(i1);
    if(threei2 < 3*i1)
    {
        int threeind1_i2 = 3*(ind1_i2+1);
        elems[threei2][threeind1_i2] = elems[threei2][threeind1_i2]+ local.GetElem(6,3);
        elems[threei2+1][threeind1_i2] = elems[threei2+1][threeind1_i2]+ local.GetElem(7,3);
        elems[threei2+2][threeind1_i2] = elems[threei2+2][threeind1_i2]+ local.GetElem(8,3);
        elems[threei2][threeind1_i2+1] = elems[threei2][threeind1_i2+1]+ local.GetElem(6,4);
        elems[threei2+1][threeind1_i2+1] = elems[threei2+1][threeind1_i2+1]+ local.GetElem(7,4);
        elems[threei2+2][threeind1_i2+1] = elems[threei2+2][threeind1_i2+1]+ local.GetElem(8,4);
        elems[threei2][threeind1_i2+2] = elems[threei2][threeind1_i2+2]+ local.GetElem(6,5);
        elems[threei2+1][threeind1_i2+2] = elems[threei2+1][threeind1_i2+2]+ local.GetElem(7,5);
        elems[threei2+2][threeind1_i2+2] = elems[threei2+2][threeind1_i2+2]+ local.GetElem(8,5);
    }

    int ind3_i2 = neighrs.indexOf(i3);
    if(threei2 < 3*i3)
    {
        int threeind3_i2 = 3*(ind3_i2+1);
        elems[threei2][threeind3_i2] = elems[threei2][threeind3_i2]+ local.GetElem(6,9);
        elems[threei2+1][threeind3_i2] = elems[threei2+1][threeind3_i2]+ local.GetElem(7,9);
        elems[threei2+2][threeind3_i2] = elems[threei2+2][threeind3_i2]+ local.GetElem(8,9);
        elems[threei2][threeind3_i2+1] = elems[threei2][threeind3_i2+1]+ local.GetElem(6,10);
        elems[threei2+1][threeind3_i2+1] = elems[threei2+1][threeind3_i2+1]+ local.GetElem(7,10);
        elems[threei2+2][threeind3_i2+1] = elems[threei2+2][threeind3_i2+1]+ local.GetElem(8,10);
        elems[threei2][threeind3_i2+2] = elems[threei2][threeind3_i2+2]+ local.GetElem(6,11);
        elems[threei2+1][threeind3_i2+2] = elems[threei2+1][threeind3_i2+2]+ local.GetElem(7,11);
        elems[threei2+2][threeind3_i2+2] = elems[threei2+2][threeind3_i2+2]+ local.GetElem(8,11);
    }

    int ind4_i2 = neighrs.indexOf(i4);
    if(threei2 < 3*i4)
    {
        int threeind4_i2 = 3*(ind4_i2+1);
        elems[threei2][threeind4_i2] = elems[threei2][threeind4_i2]+ local.GetElem(6,12);
        elems[threei2+1][threeind4_i2] = elems[threei2+1][threeind4_i2]+ local.GetElem(7,12);
        elems[threei2+2][threeind4_i2] = elems[threei2+2][threeind4_i2]+ local.GetElem(8,12);
        elems[threei2][threeind4_i2+1] = elems[threei2][threeind4_i2+1]+ local.GetElem(6,13);
        elems[threei2+1][threeind4_i2+1] = elems[threei2+1][threeind4_i2+1]+ local.GetElem(7,13);
        elems[threei2+2][threeind4_i2+1] = elems[threei2+2][threeind4_i2+1]+ local.GetElem(8,13);
        elems[threei2][threeind4_i2+2] = elems[threei2][threeind4_i2+2]+ local.GetElem(6,14);
        elems[threei2+1][threeind4_i2+2] = elems[threei2+1][threeind4_i2+2]+ local.GetElem(7,14);
        elems[threei2+2][threeind4_i2+2] = elems[threei2+2][threeind4_i2+2]+ local.GetElem(8,14);
    }

    int ind5_i2 = neighrs.indexOf(i5);
    if(threei2 < 3*i5)
    {
        int threeind5_i2 = 3*(ind5_i2+1);
        elems[threei2][threeind5_i2] = elems[threei2][threeind5_i2]+ local.GetElem(6,15);
        elems[threei2+1][threeind5_i2] = elems[threei2+1][threeind5_i2]+ local.GetElem(7,15);
        elems[threei2+2][threeind5_i2] = elems[threei2+2][threeind5_i2]+ local.GetElem(8,15);
        elems[threei2][threeind5_i2+1] = elems[threei2][threeind5_i2+1]+ local.GetElem(6,16);
        elems[threei2+1][threeind5_i2+1] = elems[threei2+1][threeind5_i2+1]+ local.GetElem(7,16);
        elems[threei2+2][threeind5_i2+1] = elems[threei2+2][threeind5_i2+1]+ local.GetElem(8,16);
        elems[threei2][threeind5_i2+2] = elems[threei2][threeind5_i2+2]+ local.GetElem(6,17);
        elems[threei2+1][threeind5_i2+2] = elems[threei2+1][threeind5_i2+2]+ local.GetElem(7,17);
        elems[threei2+2][threeind5_i2+2] = elems[threei2+2][threeind5_i2+2]+ local.GetElem(8,17);
    }

    ///////////// Llenar las posiciones de la matriz correspondientes al vertice i3 /////////////
    int threei3 = 3*i3;
    elems[threei3][0] = elems[threei3][0]+ local.GetElem(9,9);
    elems[threei3+1][0] = elems[threei3+1][0]+ local.GetElem(10,9);
    elems[threei3+2][0] = elems[threei3+2][0]+ local.GetElem(11,9);
    elems[threei3][1] = elems[threei3][1]+ local.GetElem(9,10);
    elems[threei3+1][1] = elems[threei3+1][1]+ local.GetElem(10,10);
    elems[threei3+2][1] = elems[threei3+2][1]+ local.GetElem(11,10);
    elems[threei3][2] = elems[threei3][2]+ local.GetElem(9,11);
    elems[threei3+1][2] = elems[threei3+1][2]+ local.GetElem(10,11);
    elems[threei3+2][2] = elems[threei3+2][2]+ local.GetElem(11,11);

    neighrs = GetElems(i3);
    int ind0_i3 = neighrs.indexOf(i0);
    if(threei3 < 3*i0)
    {
        int threeind0_i3 = 3*(ind0_i3+1);
        elems[threei3][threeind0_i3] = elems[threei3][threeind0_i3]+ local.GetElem(9,0);
        elems[threei3+1][threeind0_i3] = elems[threei3+1][threeind0_i3]+ local.GetElem(10,0);
        elems[threei3+2][threeind0_i3] = elems[threei3+2][threeind0_i3]+ local.GetElem(11,0);
        elems[threei3][threeind0_i3+1] = elems[threei3][threeind0_i3+1]+ local.GetElem(9,1);
        elems[threei3+1][threeind0_i3+1] = elems[threei3+1][threeind0_i3+1]+ local.GetElem(10,1);
        elems[threei3+2][threeind0_i3+1] = elems[threei3+2][threeind0_i3+1]+ local.GetElem(11,1);
        elems[threei3][threeind0_i3+2] = elems[threei3][threeind0_i3+2]+ local.GetElem(9,2);
        elems[threei3+1][threeind0_i3+2] = elems[threei3+1][threeind0_i3+2]+ local.GetElem(10,2);
        elems[threei3+2][threeind0_i3+2] = elems[threei3+2][threeind0_i3+2]+ local.GetElem(11,2);
    }

    int ind1_i3 = neighrs.indexOf(i1);
    if(threei3 < 3*i1)
    {
        int threeind1_i3 = 3*(ind1_i3+1);
        elems[threei3][threeind1_i3] = elems[threei3][threeind1_i3]+ local.GetElem(9,3);
        elems[threei3+1][threeind1_i3] = elems[threei3+1][threeind1_i3]+ local.GetElem(10,3);
        elems[threei3+2][threeind1_i3] = elems[threei3+2][threeind1_i3]+ local.GetElem(11,3);
        elems[threei3][threeind1_i3+1] = elems[threei3][threeind1_i3+1]+ local.GetElem(9,4);
        elems[threei3+1][threeind1_i3+1] = elems[threei3+1][threeind1_i3+1]+ local.GetElem(10,4);
        elems[threei3+2][threeind1_i3+1] = elems[threei3+2][threeind1_i3+1]+ local.GetElem(11,4);
        elems[threei3][threeind1_i3+2] = elems[threei3][threeind1_i3+2]+ local.GetElem(9,5);
        elems[threei3+1][threeind1_i3+2] = elems[threei3+1][threeind1_i3+2]+ local.GetElem(10,5);
        elems[threei3+2][threeind1_i3+2] = elems[threei3+2][threeind1_i3+2]+ local.GetElem(11,5);
    }

    int ind2_i3 = neighrs.indexOf(i2);
    if(threei3 < 3*i2)
    {
        int threeind2_i3 = 3*(ind2_i3+1);
        elems[threei3][threeind2_i3] = elems[threei3][threeind2_i3]+ local.GetElem(9,6);
        elems[threei3+1][threeind2_i3] = elems[threei3+1][threeind2_i3]+ local.GetElem(10,6);
        elems[threei3+2][threeind2_i3] = elems[threei3+2][threeind2_i3]+ local.GetElem(11,6);
        elems[threei3][threeind2_i3+1] = elems[threei3][threeind2_i3+1]+ local.GetElem(9,7);
        elems[threei3+1][threeind2_i3+1] = elems[threei3+1][threeind2_i3+1]+ local.GetElem(10,7);
        elems[threei3+2][threeind2_i3+1] = elems[threei3+2][threeind2_i3+1]+ local.GetElem(11,7);
        elems[threei3][threeind2_i3+2] = elems[threei3][threeind2_i3+2]+ local.GetElem(9,8);
        elems[threei3+1][threeind2_i3+2] = elems[threei3+1][threeind2_i3+2]+ local.GetElem(10,8);
        elems[threei3+2][threeind2_i3+2] = elems[threei3+2][threeind2_i3+2]+ local.GetElem(11,8);
    }

    int ind4_i3 = neighrs.indexOf(i4);
    if(threei3 < 3*i4)
    {
        int threeind4_i3 = 3*(ind4_i3+1);
        elems[threei3][threeind4_i3] = elems[threei3][threeind4_i3]+ local.GetElem(9,12);
        elems[threei3+1][threeind4_i3] = elems[threei3+1][threeind4_i3]+ local.GetElem(10,12);
        elems[threei3+2][threeind4_i3] = elems[threei3+2][threeind4_i3]+ local.GetElem(11,12);
        elems[threei3][threeind4_i3+1] = elems[threei3][threeind4_i3+1]+ local.GetElem(9,13);
        elems[threei3+1][threeind4_i3+1] = elems[threei3+1][threeind4_i3+1]+ local.GetElem(10,13);
        elems[threei3+2][threeind4_i3+1] = elems[threei3+2][threeind4_i3+1]+ local.GetElem(11,13);
        elems[threei3][threeind4_i3+2] = elems[threei3][threeind4_i3+2]+ local.GetElem(9,14);
        elems[threei3+1][threeind4_i3+2] = elems[threei3+1][threeind4_i3+2]+ local.GetElem(10,14);
        elems[threei3+2][threeind4_i3+2] = elems[threei3+2][threeind4_i3+2]+ local.GetElem(11,14);
    }

    int ind5_i3 = neighrs.indexOf(i5);
    if(threei3 < 3*i5)
    {
        int threeind5_i3 = 3*(ind5_i3+1);
        elems[threei3][threeind5_i3] = elems[threei3][threeind5_i3]+ local.GetElem(9,15);
        elems[threei3+1][threeind5_i3] = elems[threei3+1][threeind5_i3]+ local.GetElem(10,15);
        elems[threei3+2][threeind5_i3] = elems[threei3+2][threeind5_i3]+ local.GetElem(11,15);
        elems[threei3][threeind5_i3+1] = elems[threei3][threeind5_i3+1]+ local.GetElem(9,16);
        elems[threei3+1][threeind5_i3+1] = elems[threei3+1][threeind5_i3+1]+ local.GetElem(10,16);
        elems[threei3+2][threeind5_i3+1] = elems[threei3+2][threeind5_i3+1]+ local.GetElem(11,16);
        elems[threei3][threeind5_i3+2] = elems[threei3][threeind5_i3+2]+ local.GetElem(9,17);
        elems[threei3+1][threeind5_i3+2] = elems[threei3+1][threeind5_i3+2]+ local.GetElem(10,17);
        elems[threei3+2][threeind5_i3+2] = elems[threei3+2][threeind5_i3+2]+ local.GetElem(11,17);
    }

    ///////////// Llenar las posiciones de la matriz correspondientes al vertice i4 /////////////
    int threei4 = 3*i4;
    elems[threei4][0] = elems[threei4][0]+ local.GetElem(12,12);
    elems[threei4+1][0] = elems[threei4+1][0]+ local.GetElem(13,12);
    elems[threei4+2][0] = elems[threei4+2][0]+ local.GetElem(14,12);
    elems[threei4][1] = elems[threei4][1]+ local.GetElem(12,13);
    elems[threei4+1][1] = elems[threei4+1][1]+ local.GetElem(13,13);
    elems[threei4+2][1] = elems[threei4+2][1]+ local.GetElem(14,13);
    elems[threei4][2] = elems[threei4][2]+ local.GetElem(12,14);
    elems[threei4+1][2] = elems[threei4+1][2]+ local.GetElem(13,14);
    elems[threei4+2][2] = elems[threei4+2][2]+ local.GetElem(14,14);

    neighrs = GetElems(i4);
    int ind0_i4 = neighrs.indexOf(i0);
    if(threei4 < 3*i0)
    {
        int threeind0_i4 = 3*(ind0_i4+1);
        elems[threei4][threeind0_i4] = elems[threei4][threeind0_i4]+ local.GetElem(12,0);
        elems[threei4+1][threeind0_i4] = elems[threei4+1][threeind0_i4]+ local.GetElem(13,0);
        elems[threei4+2][threeind0_i4] = elems[threei4+2][threeind0_i4]+ local.GetElem(14,0);
        elems[threei4][threeind0_i4+1] = elems[threei4][threeind0_i4+1]+ local.GetElem(12,1);
        elems[threei4+1][threeind0_i4+1] = elems[threei4+1][threeind0_i4+1]+ local.GetElem(13,1);
        elems[threei4+2][threeind0_i4+1] = elems[threei4+2][threeind0_i4+1]+ local.GetElem(14,1);
        elems[threei4][threeind0_i4+2] = elems[threei4][threeind0_i4+2]+ local.GetElem(12,2);
        elems[threei4+1][threeind0_i4+2] = elems[threei4+1][threeind0_i4+2]+ local.GetElem(13,2);
        elems[threei4+2][threeind0_i4+2] = elems[threei4+2][threeind0_i4+2]+ local.GetElem(14,2);
    }

    int ind1_i4 = neighrs.indexOf(i1);
    if(threei4 < 3*i1)
    {
        int threeind1_i4 = 3*(ind1_i4+1);
        elems[threei4][threeind1_i4] = elems[threei4][threeind1_i4]+ local.GetElem(12,3);
        elems[threei4+1][threeind1_i4] = elems[threei4+1][threeind1_i4]+ local.GetElem(13,3);
        elems[threei4+2][threeind1_i4] = elems[threei4+2][threeind1_i4]+ local.GetElem(14,3);
        elems[threei4][threeind1_i4+1] = elems[threei4][threeind1_i4+1]+ local.GetElem(12,4);
        elems[threei4+1][threeind1_i4+1] = elems[threei4+1][threeind1_i4+1]+ local.GetElem(13,4);
        elems[threei4+2][threeind1_i4+1] = elems[threei4+2][threeind1_i4+1]+ local.GetElem(14,4);
        elems[threei4][threeind1_i4+2] = elems[threei4][threeind1_i4+2]+ local.GetElem(12,5);
        elems[threei4+1][threeind1_i4+2] = elems[threei4+1][threeind1_i4+2]+ local.GetElem(13,5);
        elems[threei4+2][threeind1_i4+2] = elems[threei4+2][threeind1_i4+2]+ local.GetElem(14,5);
    }

    int ind2_i4 = neighrs.indexOf(i2);
    if(threei4 < 3*i2)
    {
        int threeind2_i4 = 3*(ind2_i4+1);
        elems[threei4][threeind2_i4] = elems[threei4][threeind2_i4]+ local.GetElem(12,6);
        elems[threei4+1][threeind2_i4] = elems[threei4+1][threeind2_i4]+ local.GetElem(13,6);
        elems[threei4+2][threeind2_i4] = elems[threei4+2][threeind2_i4]+ local.GetElem(14,6);
        elems[threei4][threeind2_i4+1] = elems[threei4][threeind2_i4+1]+ local.GetElem(12,7);
        elems[threei4+1][threeind2_i4+1] = elems[threei4+1][threeind2_i4+1]+ local.GetElem(13,7);
        elems[threei4+2][threeind2_i4+1] = elems[threei4+2][threeind2_i4+1]+ local.GetElem(14,7);
        elems[threei4][threeind2_i4+2] = elems[threei4][threeind2_i4+2]+ local.GetElem(12,8);
        elems[threei4+1][threeind2_i4+2] = elems[threei4+1][threeind2_i4+2]+ local.GetElem(13,8);
        elems[threei4+2][threeind2_i4+2] = elems[threei4+2][threeind2_i4+2]+ local.GetElem(14,8);
    }

    int ind3_i4 = neighrs.indexOf(i3);
    if(threei4 < 3*i3)
    {
        int threeind3_i4 = 3*(ind3_i4+1);
        elems[threei4][threeind3_i4] = elems[threei4][threeind3_i4]+ local.GetElem(12,9);
        elems[threei4+1][threeind3_i4] = elems[threei4+1][threeind3_i4]+ local.GetElem(13,9);
        elems[threei4+2][threeind3_i4] = elems[threei4+2][threeind3_i4]+ local.GetElem(14,9);
        elems[threei4][threeind3_i4+1] = elems[threei4][threeind3_i4+1]+ local.GetElem(12,10);
        elems[threei4+1][threeind3_i4+1] = elems[threei4+1][threeind3_i4+1]+ local.GetElem(13,10);
        elems[threei4+2][threeind3_i4+1] = elems[threei4+2][threeind3_i4+1]+ local.GetElem(14,10);
        elems[threei4][threeind3_i4+2] = elems[threei4][threeind3_i4+2]+ local.GetElem(12,11);
        elems[threei4+1][threeind3_i4+2] = elems[threei4+1][threeind3_i4+2]+ local.GetElem(13,11);
        elems[threei4+2][threeind3_i4+2] = elems[threei4+2][threeind3_i4+2]+ local.GetElem(14,11);
    }

    int ind5_i4 = neighrs.indexOf(i5);
    if(threei4 < 3*i5)
    {
        int threeind5_i4 = 3*(ind5_i4+1);
        elems[threei4][threeind5_i4] = elems[threei4][threeind5_i4]+ local.GetElem(12,15);
        elems[threei4+1][threeind5_i4] = elems[threei4+1][threeind5_i4]+ local.GetElem(13,15);
        elems[threei4+2][threeind5_i4] = elems[threei4+2][threeind5_i4]+ local.GetElem(14,15);
        elems[threei4][threeind5_i4+1] = elems[threei4][threeind5_i4+1]+ local.GetElem(12,16);
        elems[threei4+1][threeind5_i4+1] = elems[threei4+1][threeind5_i4+1]+ local.GetElem(13,16);
        elems[threei4+2][threeind5_i4+1] = elems[threei4+2][threeind5_i4+1]+ local.GetElem(14,16);
        elems[threei4][threeind5_i4+2] = elems[threei4][threeind5_i4+2]+ local.GetElem(12,17);
        elems[threei4+1][threeind5_i4+2] = elems[threei4+1][threeind5_i4+2]+ local.GetElem(13,17);
        elems[threei4+2][threeind5_i4+2] = elems[threei4+2][threeind5_i4+2]+ local.GetElem(14,17);
    }

    ///////////// Llenar las posiciones de la matriz correspondientes al vertice i5 /////////////
    int threei5 = 3*i5;
    elems[threei5][0] = elems[threei5][0]+ local.GetElem(15,15);
    elems[threei5+1][0] = elems[threei5+1][0]+ local.GetElem(16,15);
    elems[threei5+2][0] = elems[threei5+2][0]+ local.GetElem(17,15);
    elems[threei5][1] = elems[threei5][1]+ local.GetElem(15,16);
    elems[threei5+1][1] = elems[threei5+1][1]+ local.GetElem(16,16);
    elems[threei5+2][1] = elems[threei5+2][1]+ local.GetElem(17,16);
    elems[threei5][2] = elems[threei5][2]+ local.GetElem(15,17);
    elems[threei5+1][2] = elems[threei5+1][2]+ local.GetElem(16,17);
    elems[threei5+2][2] = elems[threei5+2][2]+ local.GetElem(17,17);

    neighrs = GetElems(i5);
    int ind0_i5 = neighrs.indexOf(i0);
    if(threei5 < 3*i0)
    {
        int threeind0_i5 = 3*(ind0_i5+1);
        elems[threei5][threeind0_i5] = elems[threei5][threeind0_i5]+ local.GetElem(15,0);
        elems[threei5+1][threeind0_i5] = elems[threei5+1][threeind0_i5]+ local.GetElem(16,0);
        elems[threei5+2][threeind0_i5] = elems[threei5+2][threeind0_i5]+ local.GetElem(17,0);
        elems[threei5][threeind0_i5+1] = elems[threei5][threeind0_i5+1]+ local.GetElem(15,1);
        elems[threei5+1][threeind0_i5+1] = elems[threei5+1][threeind0_i5+1]+ local.GetElem(16,1);
        elems[threei5+2][threeind0_i5+1] = elems[threei5+2][threeind0_i5+1]+ local.GetElem(17,1);
        elems[threei5][threeind0_i5+2] = elems[threei5][threeind0_i5+2]+ local.GetElem(15,2);
        elems[threei5+1][threeind0_i5+2] = elems[threei5+1][threeind0_i5+2]+ local.GetElem(16,2);
        elems[threei5+2][threeind0_i5+2] = elems[threei5+2][threeind0_i5+2]+ local.GetElem(17,2);
    }

    int ind1_i5 = neighrs.indexOf(i1);
    if(threei5 < 3*i1)
    {
        int threeind1_i5 = 3*(ind1_i5+1);
        elems[threei5][threeind1_i5] = elems[threei5][threeind1_i5]+ local.GetElem(15,3);
        elems[threei5+1][threeind1_i5] = elems[threei5+1][threeind1_i5]+ local.GetElem(16,3);
        elems[threei5+2][threeind1_i5] = elems[threei5+2][threeind1_i5]+ local.GetElem(17,3);
        elems[threei5][threeind1_i5+1] = elems[threei5][threeind1_i5+1]+ local.GetElem(15,4);
        elems[threei5+1][threeind1_i5+1] = elems[threei5+1][threeind1_i5+1]+ local.GetElem(16,4);
        elems[threei5+2][threeind1_i5+1] = elems[threei5+2][threeind1_i5+1]+ local.GetElem(17,4);
        elems[threei5][threeind1_i5+2] = elems[threei5][threeind1_i5+2]+ local.GetElem(15,5);
        elems[threei5+1][threeind1_i5+2] = elems[threei5+1][threeind1_i5+2]+ local.GetElem(16,5);
        elems[threei5+2][threeind1_i5+2] = elems[threei5+2][threeind1_i5+2]+ local.GetElem(17,5);
    }

    int ind2_i5 = neighrs.indexOf(i2);
    if(threei5 < 3*i2)
    {
        int threeind2_i5 = 3*(ind2_i5+1);
        elems[threei5][threeind2_i5] = elems[threei5][threeind2_i5]+ local.GetElem(15,6);
        elems[threei5+1][threeind2_i5] = elems[threei5+1][threeind2_i5]+ local.GetElem(16,6);
        elems[threei5+2][threeind2_i5] = elems[threei5+2][threeind2_i5]+ local.GetElem(17,6);
        elems[threei5][threeind2_i5+1] = elems[threei5][threeind2_i5+1]+ local.GetElem(15,7);
        elems[threei5+1][threeind2_i5+1] = elems[threei5+1][threeind2_i5+1]+ local.GetElem(16,7);
        elems[threei5+2][threeind2_i5+1] = elems[threei5+2][threeind2_i5+1]+ local.GetElem(17,7);
        elems[threei5][threeind2_i5+2] = elems[threei5][threeind2_i5+2]+ local.GetElem(15,8);
        elems[threei5+1][threeind2_i5+2] = elems[threei5+1][threeind2_i5+2]+ local.GetElem(16,8);
        elems[threei5+2][threeind2_i5+2] = elems[threei5+2][threeind2_i5+2]+ local.GetElem(17,8);
    }

    int ind3_i5 = neighrs.indexOf(i3);
    if(threei5 < 3*i3)
    {
        int threeind3_i5 = 3*(ind3_i5+1);
        elems[threei5][threeind3_i5] = elems[threei5][threeind3_i5]+ local.GetElem(15,9);
        elems[threei5+1][threeind3_i5] = elems[threei5+1][threeind3_i5]+ local.GetElem(16,9);
        elems[threei5+2][threeind3_i5] = elems[threei5+2][threeind3_i5]+ local.GetElem(17,9);
        elems[threei5][threeind3_i5+1] = elems[threei5][threeind3_i5+1]+ local.GetElem(15,10);
        elems[threei5+1][threeind3_i5+1] = elems[threei5+1][threeind3_i5+1]+ local.GetElem(16,10);
        elems[threei5+2][threeind3_i5+1] = elems[threei5+2][threeind3_i5+1]+ local.GetElem(17,10);
        elems[threei5][threeind3_i5+2] = elems[threei5][threeind3_i5+2]+ local.GetElem(15,11);
        elems[threei5+1][threeind3_i5+2] = elems[threei5+1][threeind3_i5+2]+ local.GetElem(16,11);
        elems[threei5+2][threeind3_i5+2] = elems[threei5+2][threeind3_i5+2]+ local.GetElem(17,11);
    }

    int ind4_i5 = neighrs.indexOf(i4);
    if(threei5 < 3*i4)
    {
        int threeind4_i5 = 3*(ind4_i5+1);
        elems[threei5][threeind4_i5] = elems[threei5][threeind4_i5]+ local.GetElem(15,12);
        elems[threei5+1][threeind4_i5] = elems[threei5+1][threeind4_i5]+ local.GetElem(16,12);
        elems[threei5+2][threeind4_i5] = elems[threei5+2][threeind4_i5]+ local.GetElem(17,12);
        elems[threei5][threeind4_i5+1] = elems[threei5][threeind4_i5+1]+ local.GetElem(15,13);
        elems[threei5+1][threeind4_i5+1] = elems[threei5+1][threeind4_i5+1]+ local.GetElem(16,13);
        elems[threei5+2][threeind4_i5+1] = elems[threei5+2][threeind4_i5+1]+ local.GetElem(17,13);
        elems[threei5][threeind4_i5+2] = elems[threei5][threeind4_i5+2]+ local.GetElem(15,14);
        elems[threei5+1][threeind4_i5+2] = elems[threei5+1][threeind4_i5+2]+ local.GetElem(16,14);
        elems[threei5+2][threeind4_i5+2] = elems[threei5+2][threeind4_i5+2]+ local.GetElem(17,14);
    }
}
void CMatrix::WriteMatrix(QString path)
{
    FILE *file = fopen(path.toLatin1().constData(),"w");
    fprintf(file,"%i %i (rows cols)\n", rows, columns);
    for(int i=0;i<rows;i++)
    {
        for(int j=0;j<columns;j++)
        {
            if(elems[i][j] >= 0)
            {
                if(elems[i][j]< 10)
                    fprintf(file, "%20.15f ",elems[i][j]);
                else fprintf(file, "%20.14f ",elems[i][j]);
            }
            else
            {
                if(elems[i][j] > -10)
                    fprintf(file, "%20.14f ",elems[i][j]);
                else fprintf(file, "%20.13f ",elems[i][j]);
            }
        }
        fprintf(file,"\n");
    }
    fclose(file);
    delete file;
}
void CMatrix::WriteMatrixAndNeighbors(QString pathMatrix, QString pathNeigh)
{
    FILE *file = fopen(pathMatrix.toLatin1().constData(),"w");
    fprintf(file,"%i %i (rows cols)\n", rows, columns);
    for(int i=0;i<rows;i++)
    {
        for(int j=0;j<columns;j++)
        {
            if(elems[i][j] >= 0)
            {
                if(elems[i][j]< 10)
                    fprintf(file, "%20.15f ",elems[i][j]);
                else fprintf(file, "%20.14f ",elems[i][j]);
            }
            else
            {
                if(elems[i][j] > -10)
                    fprintf(file, "%20.14f ",elems[i][j]);
                else fprintf(file, "%20.13f ",elems[i][j]);
            }
        }
        fprintf(file,"\n");
    }
    fclose(file);

    file = fopen(pathNeigh.toLatin1().constData(),"w");
    fprintf(file,"%i %i (nRows nCols)\n", nrows, ncols);
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            if(neighbors[j*nrows+i] != -1)
                fprintf(file, "%6i ", neighbors[j*nrows+i]+1);
            else fprintf(file, "%6i ", -1);
        }
        fprintf(file,"\n");
    }
    fclose(file);
    delete file;
}


FComplexMatrix::FComplexMatrix(QObject* parent):QObject(parent)
{
    rows = cols = 0;
}
FComplexMatrix::FComplexMatrix(const FComplexMatrix &mat, QObject *parent){
    rows = mat.rows;
    cols = mat.cols;
    long n = rows*cols;
    if(cols > 0 && rows > 0)
        delete[] values;
    values = new complex<double>[n];
    for(int i = 0; i < n; i++)
        values[i] = mat.values[i];
}

FComplexMatrix::FComplexMatrix(int prows, int pcols,QObject* parent):QObject(parent)
{
    rows = prows;
    cols = pcols;
    values = new complex<double>[rows*cols];
}
FComplexMatrix::FComplexMatrix(complex<double> *pvalues, int prows, int pcols,QObject* parent):QObject(parent)
{
    rows = prows;
    cols = pcols;
    values = pvalues;
}

void FComplexMatrix::Clear()
{
    for(int i = 0; i < rows*cols; i++)
    {
        values[i].real(0);
        values[i].imag(0);
    }
}
void FComplexMatrix::Hermitize()
{
    for(int i = 0; i < rows; i++)
    {
        for(int j = i+1; j < cols; j++)
        {
            complex<double> herm = GetElem(i,j);
            herm.imag(-herm.imag());
            SetElem(j,i,herm);
        }
    }
}
complex<double> * FComplexMatrix::GetColumn(int ith)
{
    complex<double> *ithColumn = new complex<double>[rows];
    for(int i = 0; i < rows; i++)
    {
        complex<double> ithElem;
        ithElem.real(GetElem(i, ith).real());
        ithElem.imag(GetElem(i, ith).imag());

        ithColumn[i] = ithElem;
    }
    return ithColumn;
}
complex<double> FComplexMatrix::GetElem(int fila, int columna)
{
    return values[fila+columna*rows];
}
void FComplexMatrix::SetElem(int fila, int columna, complex<double> elem)
{
    values[fila+columna*rows] = elem;
}
void FComplexMatrix::WriteRealComplexMatrix(QString path)
{
    FILE *file = fopen(path.toLatin1().constData(),"w");
    fprintf(file,"%i %i (rows cols)\n", rows, cols);
    for(int i=0;i<rows;i++)
    {
        for(int j=0;j<cols;j++)
        {
            if(GetElem(i,j).real()>=0)
            {
                if(GetElem(i,j).real()<10)
                    fprintf(file, "%20.15f ",GetElem(i,j).real());
                else fprintf(file, "%20.14f ",GetElem(i,j).real());
            }
            else
            {
                if(GetElem(i,j).real()>-10)
                    fprintf(file, "%20.14f ",GetElem(i,j).real());
                else fprintf(file, "%20.13f ",GetElem(i,j).real());
            }
        }
        fprintf(file,"\n");
    }
    fclose(file);
}
void FComplexMatrix::WriteImagComplexMatrix(QString path)
{
    FILE *file = fopen(path.toLatin1().constData(),"w");
    fprintf(file,"%i %i (rows cols)\n", rows, cols);
    for(int i=0;i<rows;i++)
    {
        for(int j=0;j<cols;j++)
        {
            if(GetElem(i,j).imag()>=0)
            {
                if(GetElem(i,j).imag()<10)
                    fprintf(file, "%20.15f ",GetElem(i,j).imag());
                else fprintf(file, "%20.14f ",GetElem(i,j).imag());
            }
            else
            {
                if(GetElem(i,j).imag()>-10)
                    fprintf(file, "%20.14f ",GetElem(i,j).imag());
                else fprintf(file, "%20.13f ",GetElem(i,j).imag());
            }
        }
        fprintf(file,"\n");
    }
    fclose(file);
}


CComplexMatrix::CComplexMatrix(int prows, int pcolumns, QList<int> pneighbors,QObject* parent):QObject(parent)
{
    rows = prows;
    cols = pcolumns;
    nrows = rows/3;
    ncols = cols/3-1;
    neighbors = pneighbors;

    values = new complex<double>*[rows];
    for(int r = 0; r < rows; r++)
        values[r] = new complex<double>[cols];
}
void CComplexMatrix::Clear()
{
    for(int r = 0; r < rows; r++)
    {
        for(int c = 0; c < cols; c++)
        {
            values[r][c].real(0);
            values[r][c].imag(0);
        }
    }
}
complex<double> * CComplexMatrix::ToArray()
{
    complex<double>* result = new complex<double>[rows*cols];
    int index = 0;
    for(int i = 0 ; i < cols; i++)
    {
        for(int j = 0; j < rows; j ++)
        {
            result[index].real(values[j][i].real());
            result[index++].imag(values[j][i].imag());
        }
    }
    return result;
}
void CComplexMatrix::WriteRealComplexMatrix(QString path)
{
    FILE *file = fopen(path.toLatin1().constData(),"w");
    fprintf(file,"%i %i (rows cols)\n", rows, cols);
    for(int i=0;i<rows;i++)
    {
        for(int j=0;j<cols;j++)
        {
            if(values[i][j].real()>=0)
            {
                if(values[i][j].real()<10)
                    fprintf(file, "%20.15f ",values[i][j].real());
                else fprintf(file, "%20.14f ",values[i][j].real());
            }
            else
            {
                if(values[i][j].real()>-10)
                    fprintf(file, "%20.14f ",values[i][j].real());
                else fprintf(file, "%20.13f ",values[i][j].real());
            }
        }
        fprintf(file,"\n");
    }
    fclose(file);
    delete file;
}
void CComplexMatrix::WriteImagComplexMatrix(QString path)
{
    FILE *file = fopen(path.toLatin1().constData(),"w");
    fprintf(file,"%i %i (rows cols)\n", rows, cols);
    for(int i=0;i<rows;i++)
    {
        for(int j=0;j<cols;j++)
        {
            if(values[i][j].imag()>=0)
            {
                if(values[i][j].imag()<10)
                    fprintf(file, "%20.15f ",values[i][j].imag());
                else fprintf(file, "%20.14f ",values[i][j].imag());
            }
            else
            {
                if(values[i][j].imag()>-10)
                    fprintf(file, "%20.14f ",values[i][j].imag());
                else fprintf(file, "%20.13f ",values[i][j].imag());
            }
        }
        fprintf(file,"\n");
    }
    fclose(file);
    delete file;
}


int Sign(double x)
{
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}
QList<int> Sort(QList<int>elems)
{
    for(int i = 0; i < elems.count(); i++)
    {
        for(int j = 0; j < elems.count(); j++)
        {
            if(elems[j]>elems[i])
            {
                int temp = elems[i];
                elems[i] = elems[j];
                elems[j] = temp;
            }
        }
    }
    return elems;
}
