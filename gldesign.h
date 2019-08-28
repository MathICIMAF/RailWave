#ifndef GLDESIGN_H
#define GLDESIGN_H

#include <QGLWidget>
#include <cmath>
#include "RailNeeds.h"
#include "Surfaces.h"
#include <QWheelEvent>
#include <QKeyEvent>
#include <time.h>

using namespace std;

class GLDesign : public QGLWidget
{
    Q_OBJECT
public:
    explicit GLDesign(QWidget* qw = 0);
    GLDesign();
    void drawpolygon(RailPolygon* pol);
    void drawmesh(Mesh* mesh,RailDesigned* rail);
    void refine();
protected:
    void paintGL();
    virtual void initializeGL();
    void resizeGL(int w, int h);
    void wheelEvent(QWheelEvent *wheelEvent);
    void keyPressEvent(QKeyEvent *);
private:
    //QPointF a,b,c,d,e,f,g,h,bp,cp,dp,fp,gp,hp;
    RailDesigned *rail;
    Mesh *mesh;
    bool bmesh;
    QList<QPointF> polyg,CInitial;
    double radius,ratio;
    bool bcontrol;
    double obj_pos[3];
    double center[3];
    double scale;
    void ZoomIn();
    void ZoomOut();
    void Center();
    QPointF newMovementPosition;
    QPointF oldMovementPosition;
    void MoveDown();
    void MoveLeft();
    void MoveRight();
    void MoveUp();
    void draw_tri(double x0, double y0, double z0,
                   double x1, double y1, double z1,
                   double x2, double y2, double z2);
    int level;
};

#endif // GLDESIGN_H
