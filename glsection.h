#ifndef GLSECTION_H
#define GLSECTION_H

#include <QGLWidget>
#include <limits>
#include <math.h>
#include "Surfaces.h"
#include <QPointer>

#define WINDOW_WIDTH  250
#define WINDOW_HEIGHT 250


using namespace std;

class GLSection : public QGLWidget
{
    Q_OBJECT
public:
    explicit GLSection(QWidget *qw = 0);
    GLSection();
    void drawBorder(QList<QPointF> pts/*QList<PointS> pts*/);
    void drawMesh(Mesh* mesh);
    void drawControl(QList<PointS> pts);
    void GetTrianglesBoundaryRelation(QList<double> &x,QList<double> &y,QList<QList<int> > &triangles, QList<QList<int> > &boundaryIndices, QList<QList<GVector2> > &boundary);
    void CleanControl();
protected:
    void paintGL();
    virtual void initializeGL();
    void resizeGL(int w, int h);
    //void resizeGL(int width, int height);
private:
    QPointer<Mesh> mesh;
    bool border,bcontrol;
    bool dmesh;
    QList<PointS> control;
    QList<QPointF> frontier;
    double center[3];
    double obj_pos[3];
    void Center();
    void LocateCenter();
    void GetRadius();
    float scale;
    float radius;
    void displayMesh();
    void drawTriangle(float x0, float y0, float x1, float y1, float x2, float y2);

signals:

public slots:

};

#endif // GLSECTION_H
