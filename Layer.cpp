#include "Layer.h"

Material::Material(QObject*parent):QObject(parent)
{
    name = "Steel";//defaultMtrName;
    density = defaultDensity;
    shearVel = defaultShearVelocity;
    longVel = defaultLongitudinalVelocity;
    isNull = true;
}
Material::Material(const Material &mat,QObject*parent):QObject(parent){
    name = mat.GetName();
    density = mat.GetDensity();
    shearVel = mat.GetShearVelocity();
    longVel = mat.GetLongitudinalVelocity();
    isNull = mat.IsNull();
}

Material::Material(QString pname, double pdensity, double pshear, double plongVel,QObject*parent):QObject(parent)
{
    name = pname;
    density = pdensity;
    shearVel = pshear;
    longVel = plongVel;
    isNull = false;
}
bool Material::IsNull()const{ return isNull;}
QString Material::GetName()const{return name;}
double Material::GetDensity()const{ return density;}
void Material::SetNull(bool val){isNull = val;}
void Material::SetName(QString val){name = val; isNull = false;}
double Material::GetShearVelocity()const{ return shearVel;}
void Material::SetDensity(double val){ density = val; isNull = false;}
double Material::GetLongitudinalVelocity()const{ return longVel;}
void Material::SetShearVelocity(double val){ shearVel = val; isNull = false;}
void Material::SetLongitudinalVelocity(double val){ longVel = val; isNull = false;}
Material& Material::operator =(const Material &mat){
    this->name = mat.GetName();
    this->shearVel = mat.GetShearVelocity();
    this->density = mat.GetDensity();
    this->isNull = mat.IsNull();
    this->longVel = mat.GetLongitudinalVelocity();
    //this->setParent();
    return *this;
}

Layer::Layer(const Layer &layer,QObject*parent):QObject(parent)
{
    Material mat = layer.mtr;
    mtr = Material(mat.GetName(), mat.GetDensity(),mat.GetShearVelocity(), mat.GetLongitudinalVelocity());
    ShapeType shp = layer.shape;
    shape = shp;
}
Layer::Layer(ShapeType pShape, Material pmtr, QObject *parent):QObject(parent)
{
    mtr = pmtr;
    shape = pShape;
}
ShapeType Layer::GetShape(){return shape;}
Material Layer::GetMaterial(){return mtr;}
void Layer::SetShape(ShapeType type){shape = type;}
void Layer::SetMaterial(Material val){mtr = val; mtr.SetNull(false);}

//RailLayer::RailLayer(RailLayer &layer):Layer(RailEnum, layer.GetMaterial())
//{
//    aproxSpline = layer.GetAproxSpline();
//}
RailLayer::RailLayer(QList<PointS> pAproxSpline, Material pmtr,QObject*parent):Layer(RailEnum,pmtr,parent){aproxSpline = pAproxSpline;}
QList<PointS> RailLayer::GetAproxSpline(){return aproxSpline;}
void RailLayer::SetAproxSpline(QList<PointS> aprox){aproxSpline = aprox;}
