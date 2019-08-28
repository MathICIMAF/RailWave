#ifndef GLRAIL3D_H
#define GLRAIL3D_H

#include <QGLWidget>
#include "gvector3.h"
#include <QTimer>
#include <time.h>
#include "glsection.h"
#include "WavesUtils.h"
#include <QList>
#include <QTimer>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QGraphicsPixmapItem>
#include <QtOpenGL/QGLWidget>
#include <QWheelEvent>
#include <QKeyEvent>

using namespace std;

class GLRail3d : public QGLWidget
{
    Q_OBJECT
public:
    explicit GLRail3d(QWidget *qw);
    GLRail3d();
    void ShowShape(Mesh*sec, int sections, double spacing);
    void ZoomIn();
    void ZoomOut();
    void SetSpacing(double value);
    void SetSections(int value);
    void CleanControl();
    QList<QList<int> > triangles;
    QList<QList<int> > boundaryIndices;
    QList<QList<GVector2> > boundary;
    QList<double>x,y;
    void Animate(QList<QList<AvalAvect> > disp,QList<QList<AvalAvect> > pfDisplacements,QList<QList<AvalAvect> > plDisplacements,QList<QList<AvalAvect> > ptDisplacements);
protected:
    virtual void initializeGL();
    void paintGL();    
    void resizeGL(int width, int height);
    void wheelEvent(QWheelEvent *wheelEvent);
    void mouseMoveEvent(QMouseEvent *mouseEvent);
    void mousePressEvent(QMouseEvent *mousePressed);
    void keyPressEvent(QKeyEvent *);
private:
    QLabel* image;
    QList<QList<QList<GVector3*> > >InitDisplacements(int boundaries, int count,int sections);
    MeshCuad *mesh;
    Mesh* initial;
    float scale;
    float azimuth_,elevation_;
    void Center();
    float center[3];
    float radius,ratio;
    float obj_pos[3];
    void DisplayMesh();
    void LocateCenter();
    bool dmesh,animate,leftButton;
    GLSection *sec;
    int sections,timeDiscretization;
    double spacing,epsilon;
    QList<QList<QList<GVector3*> > > disp;
    void SelectionChanged();
    int row, col, type;
    QList<QList<AvalAvect> > displacements, fDisplacements, lDisplacements, tDisplacements;
    QTimer *timer;
    QList<QList<GVector3> > BoundaryDisplacements(QList<GVector3> allDisplacements,
                                                  QList<QList<int> > boundaryIndexes);
    void GetTrianglesBoundaryRelation();
    QPointF newMovementPosition;
    QPointF oldMovementPosition;
    void MoveDown();
    void MoveLeft();
    void MoveRight();
    void MoveUp();
signals:

private slots:
    void UpdateRowCol(int row, int col, int type);
    void timeElapsed();

};

#endif // GLRAIL3D_H
