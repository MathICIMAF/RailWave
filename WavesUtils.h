#ifndef WAVESUTILS_H
#define WAVESUTILS_H

#include <Matrix.h>
#include <gvector2.h>
#include <gvector3.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <QProgressDialog>
#include <QThread>
#include <QObject>
#include "Constants.h"
#include <complex>

using namespace std;

class AvalAvect:public QObject
{
private:
    int index;
    double autoval;
    double waveNumber;
    QList<complex<double> > autovect;

    double AVectRealMax();

public:
    AvalAvect(QObject* parent=0):QObject(parent){}
    AvalAvect(const AvalAvect& val, QObject* parent=0);
    AvalAvect(double WaveNumber, int index, double aVal, QList<complex<double> > aVect,QObject* parent=0);
    AvalAvect& operator=(const AvalAvect& val);
    int GetIndex() const;
    double GetAutoVal()const;
    double GetWaveNumber()const;
    void Normalize2AVector();
    void Normalize3AVector();
    void SetAutoVal(double aVal);
    QList<complex<double> > GetAutoVect()const;
    double AvectorNorm(complex<double> *x);
    complex<double> * ConvertAvectsToPointer();
    void AddValueToAutoVect(complex<double> value);
    void SetAutoVect(QList<complex<double> > aVect);
    complex<double> AvectorNorm(int initPos, int inc);
    void SetWaveNumber_Index(double pwaveNumber, int pindex);
    void ChangeValueInAutoVect(complex<double> value, int pos);
    QList<GVector3> GetDisplacements(double z, double epsilon, int t);
};

class CVertex:public QObject
{
private:
    int sectionId;
    double x, y, z;
    GVector3 *displacement;

public:
    CVertex(double pX, double pY, double pZ, int pSectionId,QObject* parent=0);
    CVertex(QObject* parent=0):QObject(parent){}
    CVertex(const CVertex &vert,QObject* parent=0);
    ~CVertex(){}
    CVertex& operator=(const CVertex &vert);
    double X()const;
    double Y()const;
    double Z()const;
    int GetSectionId()const;
    void X(double newX);
    void Y(double newY);
    void Z(double newZ);
    GVector3 * Displacement()const;
    void Displacement(GVector3 *value);
};

class Cuad:public QObject
{
private:
    CVertex vertices[4];

public:
    Cuad(CVertex v0, CVertex v1, CVertex v2, CVertex v3,QObject* parent=0);
    Cuad(const Cuad& cuad,QObject* parent=0);
    Cuad& operator=(const Cuad& cuad);
    CVertex Vertice(int index)const;
    void Vertice(int index, CVertex newValue);
};

class Section:public QObject
{
private:
    int id;
    QList<CVertex *> vertices;

public:
    Section(QList<CVertex *> pVertices, int id,QObject* parent=0);

    ~Section()
    {
        for(int i = 0; i < vertices.count(); i++)
            delete vertices.at(i);
    }

    int VerticesCount();
    CVertex * Vertice(int index);
};

class MeshCuad:public QObject
{
private:
    double deltaZ;
    int numberSections;
    QList<QList<Section *> > sections;

    void AdjustDeltha();
    void AdjustSections();
    void CreateMesh(QList<QList<GVector2> > pbound, QList<QList<QList<GVector3 *> > > pdisplacements);

public:
    QList<QList<int> > indexes;

    MeshCuad(QList<QList<GVector2> > bound, QList<QList<int> > boundaryIndexes, double deltaZ,
             int numberSections, QList<QList<QList<GVector3 *> > > displacements,QObject*parent=0);

    ~MeshCuad()
    {
        for(int i = 0; i < sections.count(); i++)
        {
            for(int j = 0; j < sections[i].count(); j++)
                delete sections[i][j];
        }
    }

    int SectionsCount();
    int BoundariesCount();
    void NoDisplacements();
    QList<Cuad> CrossSection();
    void SectionsNumber(int val);
    void DelthaValue(double val);
    double DelthaZ(){return deltaZ;}
    Section * BoundaryISectionJ(int i, int j);
    QList<Cuad> BoundaryICuadsSectionJ(int i, int j);
    void SetBoundaryISectionJ(int i, int j, Section *value);
    void ChangeDisplacements(QList<QList<QList<GVector3> > > newDisplacements);
};

class CuadraticVertexData
{
public:
    double x,y;
    int index;

    CuadraticVertexData(){}
    CuadraticVertexData(int pindex, double px, double py)
    {
        x = px;
        y = py;
        index = pindex;
    }
};

class WorkerThread : public QThread
{
    Q_OBJECT

    public:
       bool symmetry;
       double secTime;
       int autovalNumber;
       QList<double> NXIh;
       FMatrix vf, fM, vg;
       ProblemType problemType;
       QList<int> illConditionatedRowsChis;
       QList<double> x, symmetricalX, y, symmetricalY;
       FMatrix lVf, lF, lVg, fVf, fF, fVg, tVf, tF, tVg;
       QList<QList<int> > triangles, symmetricalTriangles;
       QList<int> xZeroBoundaryIndices, yZeroBoundaryIndices;
       double longVel, shearVel, density, step, maxWavesNumber;
       QList<QList<AvalAvect> > eigenValsVects, longitudinalEigens, torsionalEigens, flexionalEigens;

       WorkerThread(){}

       ~WorkerThread(){}

    private:
       double young, poisson;

       void doBeattieWork();
       void doQuadraticWork();
       void doSymmetricWork();
       void doLinealCompactedWork();
       void doLinealCGCompactedWork();
       void doQuadraticCompactedWork();
       void doQuadraticCGCompactedWork();

       int IsReplica(int index);
       QList<int> PostProcessResults();
       int Contains(QList<int> elems, int elem);
       void ShiftOneLeft(int pos, int steps, int size);
       int Contains(QList<QList<int> > edges, int v0, int v1);
       int * ConvertToFortranParam(QList<int> neighbors,int rows, int cols);
       void PostProcessEigenValVects(QList<int> deleted, int steps, int size);
       QList<int> GetNeighbors(QList<QList<int> > triangles, int n, int &ncols);
       FMatrix Diff(int rows, int cols, QList<QList<AvalAvect> > eigensValsVects);
       QList<double> XTimesGe(FMatrix *m, int row, QList<int> elem, int size, bool transpose);
       void SortRealNorm(complex<double> *w, QList<QList<AvalAvect> > &eigenvals,
                         int pos, complex<double> *vr, int maxNumber, int n, double waveNumber);
       void SortRealNorm(complex<double> *w, QList<QList<AvalAvect> > &eigenvals, int pos,
                         complex<double> *vr, int maxNumber, int n,
                         double waveNumber, QList<int> indicesToDelete, QList<int> symmetryElems);
       QList<int> GetQuadraticNeighbors(QList<QList<int> > quadraticTriangles, int n, int &ncols);
       double MatrixK1K2K3MT3N(double coord[], double f[], FMatrix *k1, FMatrix *k2, FMatrix *k3, FMatrix *m);
       double MatrixK1K2K3MT6N(double coord[], double f[], FMatrix *k1, FMatrix *k2, FMatrix *k3, FMatrix *m);
       QList<complex<double> > XTimesGe(FComplexMatrix *m, int row, QList<int> elem, int size, bool transpose);

    protected:
       virtual void run();

    signals:
       void finished();
       void posChanged(int pos);
       void maxIteration(int size);
};

class IOThread : public QThread
{
    Q_OBJECT

    public:
        bool save;
        int error;
        QString path;
        int sectionsCount;
        QList<double> chi;
        QString materialName;
        ShapeType structureType;
        QList<QList<int> > bIndices;
        QList<QList<GVector2> >boundaries;
        FMatrix f, fF, lF, tF, vF, fvF, lvF, tvF, vg, fvg, lvg, tvg;
        double shearVel, longVel, density, delthaZ, ir, thickness, ratio;
        QList<QList<AvalAvect> > displacements, fdisplacements, ldisplacements, tdisplacements;
        int maxCoordXFvsVF, maxCoordYFvsVF, maxCoordXFvsVG, maxCoordYFvsVG, maxCoordXKvsF, maxCoordYKvsF;

        IOThread();
        IOThread(QString ppath, QString pmaterialName, double pshearVel, double plongVel, double pdensity,
                 FMatrix pf, FMatrix pfF, FMatrix plF, FMatrix ptF,
                 FMatrix pvF, FMatrix pfvF, FMatrix plvF, FMatrix ptvF,
                 FMatrix pvg, FMatrix pfvg, FMatrix plvg, FMatrix ptvg,
                 QList<QList<AvalAvect> > pdisplacements,
                 QList<QList<AvalAvect> > pfdisplacements,
                 QList<QList<AvalAvect> > pldisplacements,
                 QList<QList<AvalAvect> > ptdisplacements,
                 int pmaxCoordXFvsVF, int pmaxCoordYFvsVF,
                 int pmaxCoordXFvsVG, int pmaxCoordYFvsVG,
                 int pmaxCoordXKvsF, int pmaxCoordYKvsF,
                 double delthaZ, int psectionsCount,
                 QList<QList<GVector2> > pboundaries, QList<QList<int> > pbIndices,
                 ShapeType pstructureType, double pir, double pthickness, double pratio);

       ~IOThread(){}

    private:
       void SaveExample();
       void LoadExample();
       void SaveToFile(QTextStream &stream, FMatrix f, FMatrix vF, FMatrix vg,
                       QList<QList<AvalAvect> > displacements);
       int ReadFromFile(QTextStream &out, FMatrix &f, FMatrix &vF, FMatrix &vg,
                         QList<QList<AvalAvect> > &displacements, QList<double> &chi);

    protected:
       virtual void run();

    signals:
       void finished();
       void posChanged(int pos);
       void maxIteration(int size);
};

double Abs(double val);
QList<QList<int> > ReadNeighbors(QString path);
void WriteVector(QString path, int *x, int size);
void WriteVector(QString path, QList<int> elems);
QList<double> ToDoubleList(double *S, int sCount);
void WriteVector(QString path, QList<double> elems);
FMatrix DecompactMatrix(CMatrix A, QString pathToSave);
FMatrix * ReadCMatrix(int *rows, int *cols, QString path);
void SaveElems(QList<QList<int> > triangles, QString path);
double Distance(double x1, double y1, double x2, double y2);
void WriteVector(QString path, complex<double> *x, int size);
void WriteVector(QString path, int *x, int vertexNum, int cols);
void Copy(complex<double> *from, complex<double> *to, int size);
void SaveVertices(QList<double> x, QList<double> y, QString path);
QList<complex<double> > ToComplexList(complex<double> *array, int size);
FComplexMatrix* ReadCMatrix(int *rows, int *cols, QString pathR, QString pathI);
void ScaleCoordinates(int minimumDifference, QList<double> &x, QList<double> &y);
FMatrix DecompactMatrix(QString pathToDecompact, QString neighPath, QString pathToSave);
double GetMinimumDifference(QList<QList<int> > triangles, QList<double> x, QList<double> y);
FComplexMatrix DecompactCMatrix(QString pathToDecompactR, QString pathToDecompactI,
                                QString neighPath, QString pathToSaveR, QString pathToSaveI);







#endif // WAVESUTILS_H
