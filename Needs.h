#ifndef NEEDS_H
#define NEEDS_H

#include "gvector2.h"
#include "gvector3.h"
#include "A48/a48.h"
#include <QString>
#include <QStringList>
#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include "Constants.h"


using namespace std;
using namespace A48;

struct Point
{
  GVector3 g;    // 3D position

  Point() {}
  Point(GVector3 pos) : g(pos) {}
  Point(double x, double y, double z) { g[0] = x; g[1] = y; g[2] = z; }
};
class Pvert : public Vertex
{
public:
  Point a;
  Pvert(QObject*parent=0) : Vertex(parent), a(0, 0, 0) {}
  static Pvert * cast(Vertex *v) { return static_cast<Pvert *>(v); }
};

////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Declarations /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
bool Equals(Vertex *v0, Vertex *v1);
QList<Vertex *> Get1Neighbors(Vertex *v);
QString GetDirectoryAddress(QString path);
void WriteOFFMesh(Mesh *mesh, QString path);
QList<QList<Vertex *> > GetBoundary(Mesh *s);
int Contains(QList<Edge *> vect, Edge *hedge);
void WriteOterosMesh(Mesh *mesh, QString path);
int Contains(QList<Pvert *> boundary, Pvert *vert);
double GetRadius(QList<double> x, QList<double> y);
int Contains(QList<Vertex *> boundary, Vertex *vert);
QList<Vertex *> Revert(QList<Vertex *> list, int index);
QList<Vertex *> Shift(QList<Vertex *> list, int beginIndex);
QPoint Contains(QList<QList<Vertex *> > boundary, Vertex *vert);
QList<QList<Vertex *> > SortBoundaries(QList<QList<Vertex *> > list);
void WriteBoundary(QString path, QList<int> index, QList<Vertex *> boundary);
QList<Vertex *> GetDiff(QList<Vertex *> allBoundVertices, QList<Vertex *> bound);
void WriteDATMesh(Mesh *mesh, QString path, ShapeType shape, double ir, double thickness);


#endif // NEEDS_H
