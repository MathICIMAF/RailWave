#ifndef GVECTOR2_H
#define GVECTOR2_H

/* Declarations */
class GVector2;

/* Dependencies */
#include <math.h>
#include <iostream>
#include <QObject>

/* Definitions */

/** \brief Vector in 3D. 
*	\ingroup Vectors
*/
class GVector2:public QObject
{
public:
	/** \brief Create a vector v=(0,0). */
        GVector2(QObject *parent=0);

	/** \brief Create a vector v=(ax,ay). */
        GVector2(double ax, double ay,QObject* parent=0);

	/** \brief Create a vector w=v(ax,ay). */
        GVector2(const GVector2& v,QObject* parent = 0);

        GVector2& operator=(const GVector2& v);

	/**Calcula la distancia euclideana*/
	double EuclideanDistance(const GVector2& v);
	/** \brief Return the component i={0,1} (ax=v[0], ay=v[1]) of the 
	vector. */
	double& operator [] (int i)
	{
		if (i == 0) return ax;
		return ay;
	}

	/** \brief Return the component i={0,1} (ax=v[0], ay=v[1]) of the 
	vector. */
	double operator [] (int i) const
	{
		if (i == 0) return ax;
		return ay;
	}

	/** \brief Return true if the coordinates of the vectors v0 and v1 are equal. */
	bool operator == (const GVector2& v)
	{
		return (ax == v.ax && ay == v.ay );
	}

	/** \brief Return true if the vector is equal to (0,0). */
        bool is_null() const
        {
                return (ax == 0.0 && ay == 0.0);
        }

		
	/** \brief Compute the Euclidean norm of the vector. */
        double norm() const;

        /** \brief Compute the squared Euclidean norm of the vector. */
        double norm2() const;

        /** \brief Normalize a vector. Divide each component by the Euclidean norm. */
        void normalize();

        /** \brief Compute the dot product between two vectors v0 and v1. */
        static double dot(const GVector2& v0, const GVector2& v1);

        /** \brief Compute the cross product between two vectors v0 and v1. */
	static GVector2* cross(const GVector2& v1, const GVector2& v2);
	
	/** \brief Compute the cosine of the angle between to vectors v1 and v2. (v1 dot v2 = |v1||v2|cos(alpha)) */
	static double cos(const GVector2& v1, const GVector2& v2);

        /** \brief Compute the sum (v=v0+v1) of two vectors v0 and v1. */
        static GVector2* add(const GVector2& v0, const GVector2& v1) ;

        /** \brief Compute the difference (v=v1-v0) of two vectors v1 and v0. */
        static GVector2* subs(const GVector2& v1, const GVector2& v0) ;

        /** \brief Compute the product (w=v*s) of a vector v and an scalar s. */
        static GVector2* mult(const GVector2& v, double s) ;

        /** \brief Compute the quotient (w=v/s) of a vector v and an scalar s.
        *       \return w=v/s if s!=0, 0 in other case.
        */
        static GVector2* div(const GVector2& v, double s) ;

        /// First vector component.
        double ax;
        /// Second vector component.
        double ay;

		GVector2 operator + (const GVector2 &v)
		{
			return *add(*this, v);
		}
		GVector2 operator - (const GVector2 &v)
		{
			return *subs(*this, v);
		}
		GVector2 operator *(const double s)
		{
			return *mult(*this, s);
		}
		GVector2 operator /(const double s)
		{
			return *div(*this, s);
		}

};

#endif
