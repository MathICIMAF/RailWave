#ifndef SURFACES_H
#define SURFACES_H

#include "Needs.h"
#include "RailNeeds.h"

#include <A48/a48.h>
#include "gvector3.h"

using namespace std;
using namespace A48;


class Rail : public virtual Surface,public QObject
{
public:
    int *fcs;
    int np, nf;
    double *pts;
    double *pts2;
    double thresh;
    int target_level;
    BinTreePoligon treePoligon;

    Rail(QObject* parent=0):QObject(parent){}
    Rail(QList<PointS> puntos, int nPoints, int nFaces,QObject* parent=0);

    ~Rail();

    float elenght(Edge *e);
    float ref_rank(Edge *e);
    float simpl_rank(Vertex *w);
    void sample(int i, Vertex *v);
    void sample(Edge *e, Vertex *v);
    void sample(Face *f, Vertex *v);
    void base_mesh(int *np, int **fcs, int *nf);
    Vertex * new_vertex() {  return new Pvert(); }//void
    void del_vertex(Vertex *v) { delete(Pvert::cast(v)); }
};

class RailDesigned:public virtual Surface{
public:
    int *fcs;
    int np, nf;
    double *pts;
    double *pts2;
    double thresh;
    int target_level;
    BinTreePoligon treePoligon;
    RailDesigned(RailPolygon*pol, int nPoints, int nFaces);

    ~RailDesigned();

    float elenght(Edge *e);
    float ref_rank(Edge *e);
    float simpl_rank(Vertex *w);
    void sample(int i, Vertex *v);
    void sample(Edge *e, Vertex *v);
    void sample(Face *f, Vertex *v);
    void base_mesh(int *np, int **fcs, int *nf);
    Vertex * new_vertex() {  return new Pvert(); }//void
    void del_vertex(Vertex *v) { delete(Pvert::cast(v)); }
};

#endif // SURFACES_H
