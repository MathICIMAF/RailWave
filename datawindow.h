#ifndef DATAWINDOW_H
#define DATAWINDOW_H

#include <QMainWindow>
//#include "graphics.h"
#include <QFutureWatcher>
#include <QtGui>
#include <QImage>
#include <iostream>
#include <QtConcurrentMap>
#include <QFuture>
#include "Layer.h"
#include <time.h>
#include <QMessageBox>
//#include "geometry.h"
#include "Utils.h"
#include "Matrix.h"
#include "WavesUtils.h"


namespace Ui {
    class datawindow;
}

class datawindow : public QMainWindow
{
    Q_OBJECT

public:
    double secTime;
    QList<double> NXIT;
    QList<int> illRows;
    int refinementSteps;
    double ir, thickness;
    double refinementThresh;
    QList<QList<int> > boundaryIndex;
    QList<QList<GVector2> > boundary;
    QList<double> x, symmetricalX, y, symmetricalY;
    QList<QList<int> > triangles, symmetricalTriangles;
    QList<int> xZeroBoundaryIndices, yZeroBoundaryIndices;
    FMatrix f, vf, vg, fF, fVf, fVg, lF, lVf, lVg, tF, tVf, tVg;
    QList<QList<AvalAvect> > eigenValVects, longitudinalEigens, torsionalEigens, flexuralEigens;

    explicit datawindow(QWidget *parent = 0);
    datawindow(int prefinemetSteps, double prefinementThresh, double ri, double thickness, QList<double> px, QList<double> py,
               QList<QList<int> > ptriangles, QList<QList<int> > pBoundary, double pshearVel,
               double plongVel, double pdensity, double pmaxWaveNumber, double pstep,
               int pcurvesNumber, bool psymmetry, ShapeType pstructure, ProblemType problemType);

    ~datawindow();

private:
    bool symmetries;
    Ui::datawindow *ui;
    ShapeType structure;
    WorkerThread worker;

    void Compute();
    QList<QList<GVector2> >GetBoundary();
    QList<QList<QList<GVector3> > > GetDisplacements(int boundaries, int count, int sections);

private slots:
    void finish();
    void on_cancelBtn_clicked();
    void closeEvent(QCloseEvent *);
    void setProgressBarPos(int pos);

signals:
    void SetProgressBar(int);
    void Cancelled(datawindow*);
    void CurvesToShow(bool, datawindow*);
    void SetHidden(int val, int max, datawindow*);
};

#endif // DATAWINDOW_H
