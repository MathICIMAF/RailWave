#ifndef GVECTOR3_H
#define GVECTOR3_H

/* Declarations */
class GVector3;

/* Dependencies */
#include <math.h>
#include <iostream>
#include <QObject>

/* Definitions */

/** \brief Vector in 3D. 
*	\ingroup Vectors
*/
class GVector3
{
public:
	/** \brief Create a vector v=(0,0,0). */
        GVector3();

	/** \brief Create a vector v=(ax,ay,az). */
        GVector3(double ax, double ay, double az);

	/** \brief Create a vector w=v(ax,ay,az). */
        GVector3(const GVector3& v);

        //GVector3 operator=(const GVector3 v);

	/**Calcula la distancia euclideana*/
        double EuclideanDistance(const GVector3& v);

	/** \brief Return the component i={0,1,2} (ax=v[0], ay=v[1], az=v[2]) of the 
	vector. */
	double& operator [] (int i)
	{
		if (i == 0) return ax;
		else if (i == 1) return ay;
		return az;
	}

	/** \brief Return the component i={0,1,2} (ax=v[0], ay=v[1], az=v[2]) of the 
	vector. */
	double operator [] (int i) const
	{
		if (i == 0) return ax;
		else if (i == 1) return ay;
		return az;
	}

	/** \brief Return true if the coordinates of the vectors v0 and v1 are equal. */
	bool operator == (const GVector3& v)
	{
		return (ax == v.ax && ay == v.ay && az == v.az);
	}

	/** \brief Return true if the vector is equal to (0,0,0). */
        bool is_null() const
        {
                return (ax == 0.0 && ay == 0.0 && az == 0.0);
        }

		
	/** \brief Compute the Euclidean norm of the vector. */
        double norm() const;

        /** \brief Compute the squared Euclidean norm of the vector. */
        double norm2() const;

        /** \brief Normalize a vector. Divide each component by the Euclidean norm. */
        void normalize();

		

        /** \brief Compute the dot product between two vectors v0 and v1. */
        static double dot(const GVector3& v0, const GVector3& v1);

        /** \brief Compute the cross product between two vectors v0 and v1. */
	static GVector3* cross(const GVector3& v1, const GVector3& v2);
	
	/** \brief Compute the cosine of the angle between to vectors v1 and v2. (v1 dot v2 = |v1||v2|cos(alpha)) */
	static double cos(const GVector3& v1, const GVector3& v2);

        /** \brief Compute the sum (v=v0+v1) of two vectors v0 and v1. */
        static GVector3* add(const GVector3& v0, const GVector3& v1) ;

        /** \brief Compute the difference (v=v1-v0) of two vectors v1 and v0. */
        static GVector3* subs(const GVector3& v1, const GVector3& v0) ;

        /** \brief Compute the product (w=v*s) of a vector v and an scalar s. */
        static GVector3* mult(const GVector3& v, double s) ;

        /** \brief Compute the quotient (w=v/s) of a vector v and an scalar s.
        *       \return w=v/s if s!=0, 0 in other case.
        */
        static GVector3* div(const GVector3& v, double s) ;

        /// First vector component.
        double ax;
        /// Second vector component.
        double ay;
	/// Third vector component.
	double az;

		GVector3 operator + (const GVector3 &v)
		{
			return *add(*this, v);
		}
		GVector3 operator - (const GVector3 &v)
		{
			return *subs(*this, v);
		}
		GVector3 operator *(const double s)
		{
			return *mult(*this, s);
		}
		GVector3 operator /(const double s)
		{
			return *div(*this, s);
		}

};

#endif
