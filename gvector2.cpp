#include "gvector2.h"

GVector2::GVector2(QObject* parent):QObject(parent)
{
	this->ax = 0.0;
	this->ay = 0.0;
}

GVector2::GVector2(double ax, double ay, QObject *parent):QObject(parent)
{
	this->ax = ax;
	this->ay = ay;
}

GVector2::GVector2(const GVector2& v,QObject* parent):QObject(parent)
{
	this->ax = v.ax;
	this->ay = v.ay;
}

GVector2& GVector2::operator =(const GVector2& v){
    this->ax = v.ax;
    this->ay = v.ay;
    return *this;
}

double GVector2::norm() const
{
	return sqrt(ax*ax + ay*ay);
}

double GVector2::norm2() const
{
	return ax*ax + ay*ay;
}

void GVector2::normalize()
{
	double nr = sqrt(ax*ax + ay*ay);
	if (nr <= 1e-8) std::cout << "in vector normalization => The Euclidean norm is less than 1e-8..." << std::endl;
	this->ax /= nr;
	this->ay /= nr;
}

double GVector2::dot(const GVector2& v0, const GVector2& v1)
{
	return v0.ax*v1.ax + v0.ay*v1.ay;
}


GVector2* GVector2::cross(const GVector2& v0, const GVector2& v1)
{
	double ax = v0.ay - v1.ay,
	       ay = v1.ax - v0.ay;
	return new GVector2(ax, ay);
}

double GVector2::cos(const GVector2& v0, const GVector2& v1)
{
	double n0 = v0.norm(),
	       n1 = v1.norm();
	return dot(v0, v1)/(n0*n1);
}

GVector2* GVector2::add(const GVector2& v0, const GVector2& v1)
{
	return new GVector2(v0.ax+v1.ax, v0.ay+v1.ay);
}

GVector2* GVector2::subs(const GVector2& v0, const GVector2& v1)
{
	return new GVector2(v1.ax-v0.ax, v1.ay-v0.ay);
}

GVector2* GVector2::mult(const GVector2& v, double s)
{
	return new GVector2(s*v.ax, s*v.ay);
}

GVector2* GVector2::div(const GVector2& v, double s)
{
	if (s == 0.0) return 0;
	return new GVector2(v.ax/s, v.ay/s);
}

double GVector2::EuclideanDistance(const GVector2 &v)
{
	return sqrt(pow(ax-v.ax,2)+pow(ay-v.ay,2));
}
