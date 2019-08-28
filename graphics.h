#ifndef GRAPHICS_H
#define GRAPHICS_H

#include <QObject>
#include "qcustomplot.h"
#include "datawindow.h"

class graphics : public QObject
{
    Q_OBJECT
public:
    explicit graphics(QObject *parent = 0);
    graphics(QList<QString> materialInfo,FMatrix pF, FMatrix pVf,
            bool modesShape, QList<int> illRows, QCustomPlot* cPlot);
    void ShowCurves(bool showC);
private:
   int fil, col;
   QList<int> illRows;
   FMatrix X, Y;
   QString measureX, measureY;
   void GetMaxMinCoords(FMatrix X, FMatrix Y);
   double autoScaleXMin, autoScaleXMax, autoScaleYMin, autoScaleYMax, autoScaleIncX, autoScaleIncY;
   void InsertInOrder(double elem, QList<double> &list);
   void InitCustomPlot(QString title, int maxCoordX, int maxCoordY);
   QCustomPlot* customPlot;
   QString xAxisLabel, yAxisLabel;
   void InitCurvesData(FMatrix pX, FMatrix pY, QColor *color, int type);
   void ShowCurvesData(FMatrix X, FMatrix Y, QColor *color, int type);
   void ShowGraphs(FMatrix X, FMatrix Y, QColor *color, int type);
   void ChangeMeasures(QString newMeasure, double inc, QString axis);
   QList<QColor> GetColors(int);
signals:
   void SelectionChanged(int row, int col, int type);

private slots:
   void selectionChanged();
   void NewSelected(QCPAbstractPlottable *);
   void mousePress(QMouseEvent *);
   void mouseWheel(QWheelEvent *);
   void graphClicked(QCPAbstractPlottable *);
   void graphMouseMove(QMouseEvent *);

};

#endif // GRAPHICS_H
