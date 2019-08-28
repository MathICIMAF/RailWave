#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "Surfaces.h"
#include "Constants.h"
#include "datawindow.h"
#include "graphics.h"
#include "dialog.h"
#include <QPointer>
#include <designrail.h>

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
    QPointer<Rail> rail;
    QPointer<Mesh> mesh;
    Mesh *initial;
    double maxChi,step,shearVel,density,longVel;
    QString materialName;
    int curvesNo;
    ProblemType probType;
    int refinementSteps;
    double refinementThresh;
    ShapeType structureType;
    FMatrix f, vF, vG, lf, lvF, lvG, ff, fvF, fvG, tf, tvF, tvG;
    QList<QList<AvalAvect> > avalsVects, lvalsVects, tvalsVects, fvalsVects;    
    QList<double> chi;
    double maxX, maxY,radius;
    QList<int> illRows;
    double secTime;
    QList<QString> materialInfo;
    graphics *gr;
    bool railLoaded;bool graphloaded;
    bool animated;
    QList<PointS> result;
    QList<QPointF> front;
    void DrawRail(QList<PointS> res);
    RailPolygon* pol;
    RailDesigned* raild;
private slots:    
    void on_pushButton_clicked();
    void on_actionCreate_Rail_triggered();
    void on_actionAbout_triggered();
    //void on_highRefRB_clicked();
    //void on_mediumRefRB_clicked();
    void on_animatePB_clicked();
    void on_actionShow_Curves_triggered();
    void on_spacingSp_valueChanged(double );
    void on_secNumbSp_valueChanged(int );
    void on_computeButt_clicked();
    void on_actionClose_Rail_triggered();
    void on_zoomOut_clicked();
    void on_zoomIn_clicked();
    void on_actionClose_triggered();
    void on_actionOpen_Rail_triggered();
    void Cancelled(datawindow *);
    void ShowDisplayMenu(bool, datawindow *);
    void ShowRail(QList<QPointF> front,Mesh* mesh,double radius,RailDesigned* rail=0,bool refine=false);
};

#endif // MAINWINDOW_H
