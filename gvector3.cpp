#include "gvector3.h"

GVector3::GVector3()
{
	this->ax = 0.0;
	this->ay = 0.0;
	this->az = 0.0;
}

GVector3::GVector3(double ax, double ay, double az/*, QObject *parent*/)//:QObject(parent)
{
	this->ax = ax;
	this->ay = ay;
	this->az = az;
}

GVector3::GVector3(const GVector3& v/*,QObject* parent*/)//:QObject(parent)
{
	this->ax = v.ax;
	this->ay = v.ay;
	this->az = v.az;
}


double GVector3::norm() const
{
	return sqrt(ax*ax + ay*ay + az*az);
}

double GVector3::norm2() const
{
	return ax*ax + ay*ay + az*az;
}

void GVector3::normalize()
{
	double nr = sqrt(ax*ax + ay*ay + az*az);
	if (nr <= 1e-8) std::cout << "in vector normalization => The Euclidean norm is less than 1e-8..." << std::endl;
	this->ax /= nr;
	this->ay /= nr;
	this->az /= nr;
}

double GVector3::dot(const GVector3& v0, const GVector3& v1)
{
	return v0.ax*v1.ax + v0.ay*v1.ay + v0.az*v1.az;
}


GVector3* GVector3::cross(const GVector3& v0, const GVector3& v1)
{
	double ax = v0.ay*v1.az - v0.az*v1.ay,
	       ay = v1.ax*v0.az - v1.az*v0.ax,
	       az = v0.ax*v1.ay - v1.ax*v0.ay;
	return new GVector3(ax, ay, az);
}

double GVector3::cos(const GVector3& v0, const GVector3& v1)
{
	double n0 = v0.norm(),
	       n1 = v1.norm();
	return dot(v0, v1)/(n0*n1);
}

GVector3* GVector3::add(const GVector3& v0, const GVector3& v1)
{
	return new GVector3(v0.ax+v1.ax, v0.ay+v1.ay, v0.az+v1.az);
}

GVector3* GVector3::subs(const GVector3& v0, const GVector3& v1)
{
	return new GVector3(v0.ax-v1.ax, v0.ay-v1.ay, v0.az-v1.az);
}

GVector3* GVector3::mult(const GVector3& v, double s)
{
	return new GVector3(s*v.ax, s*v.ay, s*v.az);
}

GVector3* GVector3::div(const GVector3& v, double s)
{
	if (s == 0.0) return 0;
	return new GVector3(v.ax/s, v.ay/s, v.az/s);
}
double GVector3::EuclideanDistance(const GVector3 &v)
{
	return sqrt(pow(ax-v.ax,2)+pow(ay-v.ay,2)+pow(az-v.az,2));
}
