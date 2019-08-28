#ifndef RAILNEEDS_H
#define RAILNEEDS_H

#include <QString>
#include <QStringList>
#include <math.h>
#include <QList>
#include "gvector2.h"
#include "binheap.h"
#include <float.h>
#include <QPoint>
#include <cmath>
#include <vector>
#include <QFile>
#include <QString>
#include <QStringList>
#include <QTextStream>

using namespace std;

class Dictionary:public QObject{
public:
        QStringList keys;
        QList<float> values;

        Dictionary(QObject* parent=0):QObject(parent){}
        Dictionary(const Dictionary& dic,QObject* parent=0);
        ~Dictionary(){}

        void Fill(QString text);
        float GetValue( QString key);
        void Add(QString key, float value);
        bool Equals(QString key1, QString key2);                
};

class PointS:public QObject
{
    public:
        double _x;
        double _y;
        double _tan;
        bool _tanInf;
        bool _tanAsigned;

        PointS(QObject*parent=0);
        PointS(const PointS& pt, QObject* parent=0);
        PointS(double x, double y,QObject* parent=0);
        PointS(double x, double y, double tan,QObject* parent=0);
        PointS& operator=( const PointS& pt);

        ~PointS(){}

        void SetTangent(double value);
};

class Poligon:public QObject
{
public:
    int _excepcion;
    bool _isSimple;
    double _epsilon;
    double _vertexCount;
    GVector2 _endVertex;
    GVector2 _initVertex;
    QList<double> _xCoord;
    QList<double> _yCoord;

    Poligon(QObject* parent=0):QObject(parent){}
    Poligon(const Poligon& pol,QObject* parent=0);
    Poligon(QList<double> xCoord, QList<double> yCoord,QObject* parent = 0);
    Poligon& operator=(const Poligon& pol);

    ~Poligon(){}

    double PointToBaseSegmentDistance(GVector2 point);
    void MaximumDistance(QList<double> xCoord, QList<double> yCoord);
};

class BinTreePoligon:public QObject
{
public:
    double _priority;
    Poligon _nodePoligon;
    BinTreePoligon *_left;
    BinTreePoligon *_right;
    Poligon _closestNeighbor;

    BinTreePoligon(QObject* parent=0);
    BinTreePoligon(const BinTreePoligon& pol, QObject* parent=0);
    BinTreePoligon(Poligon node,QObject* parent =0);
    BinTreePoligon& operator =(const BinTreePoligon& pol);

    ~BinTreePoligon();


    double ComputePriority(GVector2 point);
    double GetPriority(){ return _priority;}
    double PointToPoligonDistance(GVector2 point);
    BinTreePoligon * GetLeftChild(){return _left;}
    BinTreePoligon * GetRightChild(){return _right;}
};

/*
        Binary Heap of Poligons.
        Each item must have a property of type 'double' with the name 'GetEpsilon()'.
        The 'GetEpsilon()' property is used to compare elements.
*/
class BinHeapPoligon : public BinHeap<BinTreePoligon*>
{
public:
        BinHeapPoligon(int max_size) : BinHeap<BinTreePoligon *>(max_size){}

        /// Add a new element. Return false if the heap is full.
        virtual int add(BinTreePoligon *x);

        /// Remove the element with the minimum key. Return this element.
        virtual BinTreePoligon * remove_min();

        /// Remove the element in position 'i'.
        virtual BinTreePoligon * remove(int i);

        /// Update the position of a node. Maybe its key changed ... and you need to update.
        void update(int i);

private:

        /// This functions may be redefined in derived classes depending on T.
        virtual bool great_than(BinTreePoligon *elema, BinTreePoligon *elemb) const
        {
                return elema->GetPriority() > elemb->GetPriority();
        }

        virtual bool less_than(BinTreePoligon *elema, BinTreePoligon *elemb) const
        {
                return elema->GetPriority() < elemb->GetPriority();
        }

        virtual bool equal_to(BinTreePoligon *elema, BinTreePoligon *elemb) const
        {
                return elema->GetPriority() == elemb->GetPriority();
        }

        /// swap elements i and j
        virtual void swap(int i, int j);
};

double toDoubleLength(QString length);
double toDoubleAngle(QString angle);
QList<float> WReader();
QList<PointS> CreatePoints(Dictionary measures);
QList<PointS> DelimitCurves(int curvenumber, QList<PointS> equis);
QList<PointS> FillSpace(int totalPoints, PointS initPoint, PointS endPoint);
void CreateTriangulation(QList<PointS> result, double *pts, double *pts2, int *fcs);
QList<QPointF> CreateTriangulation(QList<QPointF> polig, double *pts, int *fcs,QList<double> ws);
QList<PointS> AproximateSpline(QList<PointS> &curva1, QList<PointS> &curva2, QList<PointS> &curva3,
                               QList<PointS> &curva4, QList<PointS> &curva5, QList<float> ws);

////////////////// Conic Spline Section //////////////////
PointS GetQi(PointS p1, PointS p2);
//En este metodo, como se necesitan 5 puntos si pos esta entre las primeras 2 posiciones o entre
//las ultimas 2 posiciones se tomara la lista de puntos como ciclica.
double GetAutomaticTangent(QList<PointS> coord, int pos);
double GetShoulderTan(PointS initPoint, PointS endPoint);
int SubdivideEdges(QList<PointS> &result, int iterationNumber,
                   PointS initPoint, PointS endPoint, double w);
PointS IntersectionPoint(PointS p1, PointS p2, PointS p3, PointS p4);
double GetAutomaticParameter(PointS p1, PointS pq, PointS p2, PointS p);
QList<PointS> ConicSpline(int iterationNumber, QList<PointS> &coord,

                          QList<float> ws, int wIndex, int &addedPoints);
class RailPolygon: public QObject{
public:
    RailPolygon(QObject* parent = 0):QObject(parent){}
    RailPolygon(double hd,double fd,double bd,double hw,double w,
                double bw,double r1h,double r2h, double rhb, double alpha,
                double beta, double theta,QObject*parent = 0);
    void InitialPolygon();
    QList<QPointF> C(){return CInitial;}
    QList<QPointF> Polygon(){return polig;}
    double GetRadius(){return radius;}
    QPointF center;
    void SubdividePolygon(int k);
    QList<double> ws;
private:
    vector< vector<double> > puntos;
    double radius;
    double alpha, beta, theta, r1h, r2h, rhb;    
    QList<QPointF> polig,CInitial;
    QPointF a,b,c,d,e,f,g,h,bp,cp,dp,fp,gp,hp;
    double norm(QPointF p);
    QList<QPointF> Bezcircle(QPointF A, QPointF B, QPointF C, double r,double &w);
    QList<QPointF> SubdividePoints(QPointF A, QPointF B, QPointF C, int iter, double w,vector<double> &waux);


};


/*void Parameters(double hd,double fd,double bd,double hw,double w,
                    double bw,double r1h,double r2h, double rhb, double alpha,
                    double beta, double theta);*/


#endif // RAILNEEDS_H
