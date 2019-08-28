#include "Surfaces.h"

Rail::Rail(QList<PointS> puntos, int nPoints, int nFaces,QObject* parent):Surface(),QObject(parent)
{
    np = nPoints;
    nf = nFaces;
    pts = new double[np*3];
    pts2 = new double[np*2];
    fcs = new int[nf*3];
    CreateTriangulation(puntos, pts, pts2, fcs);

    target_level = 1;
    thresh = 0.008;
    QList<double> x;
    QList<double> y;

    for(int i=0;i<puntos.count();i++)
    {
        x.append(puntos.at(i)._x);
        y.append(puntos.at(i)._y);
    }
    x.append(puntos.at(0)._x);
    y.append(puntos.at(0)._y);

    treePoligon = BinTreePoligon(Poligon(x,y));
};
//devuelve la longitud de la arista
float Rail::elenght(Edge* e)
{
    Pvert* v0 = Pvert::cast(e->org());
    Pvert* v1 = Pvert::cast(e->dst());
    GVector3 p = (v1->a.g - v0->a.g);
    return p.norm();
};
//devuelve el grado de la arista de division
float Rail::ref_rank(Edge *e)
{
    Pvert* v0 = Pvert::cast(e->org());
    Pvert* v1 = Pvert::cast(e->dst());
    GVector3 p = (v0->a.g + v1->a.g) * 0.5;

    double y = treePoligon.PointToPoligonDistance(GVector2(p.ax,p.ay));

    if (qAbs(y) < thresh)
      return (e->level() < target_level + 1)? 1 : -1;
    else
      return (e->level() < target_level - 1)? 1 : -1;
};
//Devuelve el grado del vertice de refinamiento
float Rail::simpl_rank(Vertex *w)
{
    Pvert* v = Pvert::cast(w);

    double y = treePoligon.PointToPoligonDistance(GVector2(v->a.g.ax,v->a.g.ay));

    if (qAbs(y) < thresh)
    return (v->level()-1 > target_level + 1)? 1 : -1;
    else
    return (v->level()-1 > target_level - 1)? 1 : -1;
};
//Deja en v el k-esimo vertice
void Rail::sample(int k, Vertex* v)
{
    int i = k * 2;
    GVector3 p(pts2[i], pts2[i+1], 0);
    Pvert* pv = Pvert::cast(v);
    pv->a = Point(p);
};
//Deja en v el punto medio de la arista e
void Rail::sample(Edge* e, Vertex* v)
{
    Pvert* v0 = Pvert::cast(e->org());
    Pvert* v1 = Pvert::cast(e->dst());

    GVector3 p =  (v0->a.g + v1->a.g)*0.5 ;

    if(e->is_bdry())
    {
        GVector2 g(p.ax,p.ay);
        treePoligon.PointToPoligonDistance(g);
        Poligon closestPoligon = treePoligon._closestNeighbor;

        GVector2 left;
        GVector2 right;
        //Sorting vertex
        if(closestPoligon._initVertex.ax < closestPoligon._endVertex.ax)
        {
            left = closestPoligon._initVertex;
            right = closestPoligon._endVertex;
        }
        else if(closestPoligon._endVertex.ax < closestPoligon._initVertex.ax)
        {
            left = closestPoligon._endVertex;
            right = closestPoligon._initVertex;
        }
        else
        {
            if(closestPoligon._initVertex.ay < closestPoligon._endVertex.ay)
            {
                left = closestPoligon._initVertex;
                right = closestPoligon._endVertex;
            }
            else if(closestPoligon._endVertex.ay < closestPoligon._initVertex.ay)
            {
                left = closestPoligon._endVertex;
                right = closestPoligon._initVertex;
            }
            //Los dos extremos del segmentos son el mismo punto
            else
            {
                GVector3 *closest = new GVector3(closestPoligon._endVertex.ax, closestPoligon._endVertex.ay, 0);
                Pvert*pv = Pvert::cast(v);
                pv->a = Point(*closest);
            }
        }

        GVector2 qMenosP = *right.subs(left, right);
        double normaQMenosP = qMenosP.norm();

        GVector2 p0MenosP = *g.subs(left, g);

        double prodEsc = p0MenosP.dot(p0MenosP, qMenosP);
        double t = prodEsc/pow(normaQMenosP,2);

        GVector2 r = *left.add(left,*(qMenosP.mult(qMenosP,t)));

        GVector3 closest(r.ax, r.ay, 0);
        Pvert*pv = Pvert::cast(v);
        pv->a = Point(closest);
    }
    else
    {
        Pvert *pv = Pvert::cast(v);
        pv->a = Point(p);
    }
};
//Devuelve el baricentro de la cara
void Rail::sample(Face *f, Vertex* v)
{
    Pvert* v0 = Pvert::cast(f->vertex(0));
    Pvert* v1 = Pvert::cast(f->vertex(1));
    Pvert* v2 = Pvert::cast(f->vertex(2));
    GVector3 p = (v0->a.g + v1->a.g + v2->a.g) / 3.0;
    Pvert *pv = Pvert::cast(v);
    pv->a = Point(p);
};
void Rail::base_mesh(int *pnp, int **pfcs, int *pnf)
{
    *pnp = np;
    *pfcs = fcs;
    *pnf = nf;
};

Rail::~Rail(){
     delete[] pts; delete[] fcs; delete[] pts2;
}

RailDesigned::RailDesigned(RailPolygon* pol,/*QList<QPointF> C, QList<QPointF> points,*/ int nPoints, int nFaces){
    np = nPoints;
    nf = nFaces;
    pts = new double[np*3];
    fcs = new int[nf*3];    
    QList<QPointF> aditionals = CreateTriangulation(pol->C(),pts,fcs,pol->ws);
    //double xx = pts[6], yy = pts[7], zz = pts[8];
    target_level = 1;
    thresh = 0.3;
    QList<double> x;
    QList<double> y;
    QList<QPointF> points = pol->Polygon();
    for(int i=0;i<points.count();i++)
    {
        x.append(points.at(i).x());
        y.append(points.at(i).y());
    }
    /*for(int i = 0; i < aditionals.count(); i++){
        x.append(aditionals[i].x());
        y.append(aditionals[i].y());
    }*/
    treePoligon = BinTreePoligon(Poligon(x,y));
}
float RailDesigned::elenght(Edge *e){
    Pvert* v0 = Pvert::cast(e->org());
    Pvert* v1 = Pvert::cast(e->dst());
    GVector3 p = (v1->a.g - v0->a.g);
    return p.norm();
}

float RailDesigned::ref_rank(Edge *e){
    Pvert* v0 = Pvert::cast(e->org());
    Pvert* v1 = Pvert::cast(e->dst());
    GVector3 p = (v0->a.g + v1->a.g) * 0.5;

    double y = treePoligon.PointToPoligonDistance(GVector2(p.ax,p.ay));

    if (qAbs(y) < thresh)
      return (e->level() < target_level + 1)? 1 : -1;
    return (e->level() < target_level - 1)? 1 : -1;
}

float RailDesigned::simpl_rank(Vertex *w){
    Pvert* v = Pvert::cast(w);

    double y = treePoligon.PointToPoligonDistance(GVector2(v->a.g.ax,v->a.g.ay));

    if (qAbs(y) < thresh)
        return (v->level()-1 > target_level + 1)? 1 : -1;
    return (v->level()-1 > target_level - 1)? 1 : -1;
}

void RailDesigned::sample(Edge *e, Vertex *v){
    Pvert* v0 = Pvert::cast(e->org());
    Pvert* v1 = Pvert::cast(e->dst());

    GVector3 p =  (v0->a.g + v1->a.g)*0.5 ;

    if(e->is_bdry())
    {
        GVector2 g(p.ax,p.ay);
        treePoligon.PointToPoligonDistance(g);
        Poligon closestPoligon = treePoligon._closestNeighbor;

        GVector2 left;
        GVector2 right;
        //Sorting vertex
        if(closestPoligon._initVertex.ax < closestPoligon._endVertex.ax)
        {
            left = closestPoligon._initVertex;
            right = closestPoligon._endVertex;
        }
        else if(closestPoligon._endVertex.ax < closestPoligon._initVertex.ax)
        {
            left = closestPoligon._endVertex;
            right = closestPoligon._initVertex;
        }
        else
        {
            if(closestPoligon._initVertex.ay < closestPoligon._endVertex.ay)
            {
                left = closestPoligon._initVertex;
                right = closestPoligon._endVertex;
            }
            else if(closestPoligon._endVertex.ay < closestPoligon._initVertex.ay)
            {
                left = closestPoligon._endVertex;
                right = closestPoligon._initVertex;
            }
            //Los dos extremos del segmentos son el mismo punto
            else
            {
                GVector3 *closest = new GVector3(closestPoligon._endVertex.ax, closestPoligon._endVertex.ay, 0);
                Pvert*pv = Pvert::cast(v);
                pv->a = Point(*closest);
            }
        }

        GVector2 qMenosP = *right.subs(left, right);
        double normaQMenosP = qMenosP.norm();

        GVector2 p0MenosP = *g.subs(left, g);

        double prodEsc = p0MenosP.dot(p0MenosP, qMenosP);
        double t = prodEsc/pow(normaQMenosP,2);

        GVector2 r = *left.add(left,*(qMenosP.mult(qMenosP,t)));

        GVector3 closest(r.ax, r.ay, 0);
        Pvert*pv = Pvert::cast(v);
        pv->a = Point(closest);
    }
    else
    {
        Pvert *pv = Pvert::cast(v);
        pv->a = Point(p);
    }
}

void RailDesigned::sample(Face *f, Vertex *v){
    Pvert* v0 = Pvert::cast(f->vertex(0));
    Pvert* v1 = Pvert::cast(f->vertex(1));
    Pvert* v2 = Pvert::cast(f->vertex(2));
    GVector3 p = (v0->a.g + v1->a.g + v2->a.g) / 3.0;
    Pvert *pv = Pvert::cast(v);
    pv->a = Point(p);
}

void RailDesigned::base_mesh(int *pnp, int **pfcs, int *pnf){
    *pnp = np;
    *pfcs = fcs;
    *pnf = nf;
}
RailDesigned::~RailDesigned(){
    delete[] pts; delete[] fcs; delete[] pts2;
}

void RailDesigned::sample(int k, Vertex *v){
    //Provisional
    int i = k * 3;
    GVector3 p(pts[i], pts[i+1], pts[i+2]);
    Pvert* pv = Pvert::cast(v);
    pv->a = Point(p);
}
