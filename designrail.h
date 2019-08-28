#ifndef DESIGNRAIL_H
#define DESIGNRAIL_H

#include <QMainWindow>
#include <QMessageBox>
#include "RailNeeds.h"
#include "Surfaces.h"

namespace Ui {
    class DesignRail;
}

class DesignRail : public QMainWindow
{
    Q_OBJECT

public:
    explicit DesignRail(QWidget *parent = 0);
    ~DesignRail();

private:
    Ui::DesignRail *ui;
    QList<QPointF> front;
    Mesh* mesh;
    RailDesigned* rail;
    RailPolygon *pol;    
private slots:
    void on_pushButton_clicked();
    void on_actionClose_triggered();
    void on_actionSave_Rail_triggered();
    void on_refineB_clicked();
    void on_drawmesh_clicked();
    void on_drawPB_clicked();
    void on_HDlineEdit_textChanged(QString text);    
signals:
    void AddThisRail(QList<QPointF> front, Mesh* mesh,double radius,RailDesigned*,bool refine);
};

#endif // DESIGNRAIL_H
