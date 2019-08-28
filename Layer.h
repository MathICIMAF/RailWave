#ifndef LAYER_H
#define LAYER_H
#include <QString>
#include "RailNeeds.h"
#include "Constants.h"

class Material:public QObject
{
private:
    QString name;
    bool isNull;
    double density;
    double shearVel;
    double longVel;

public:
    Material(QObject*parent=0);
    Material(const Material& mat,QObject* parent=0);
    Material(QString pname, double pdensity, double pshear, double plongVel,QObject*parent=0);

    bool IsNull()const;
    QString GetName()const;
    double GetDensity()const;
    void SetNull(bool val);
    void SetName(QString val);
    double GetShearVelocity()const;
    void SetDensity(double val);
    double GetLongitudinalVelocity()const;
    void SetShearVelocity(double val);
    void SetLongitudinalVelocity(double val);
    Material& operator=(const Material&);
};

class Layer:public QObject
{
    private:
        Material mtr;
        ShapeType shape;

    public:
        Layer(QObject*parent=0):QObject(parent){}
        Layer(const Layer &,QObject* parent=0);
        Layer(ShapeType pShape, Material pmtr,QObject*parent=0);

        ShapeType GetShape();
        Material GetMaterial();
        void SetShape(ShapeType type);
        void SetMaterial(Material val);
};

class RailLayer : public Layer
{
private:
    QList<PointS> aproxSpline;

public:
    RailLayer(QObject*parent=0):Layer(parent){}
    // RailLayer(RailLayer&);
    RailLayer(QList<PointS> pAproxSpline, Material pmtr,QObject*parent=0);

    ~RailLayer(){}

    QList<PointS> GetAproxSpline();
    void SetAproxSpline(QList<PointS> aprox);
};

#endif // LAYER_H
