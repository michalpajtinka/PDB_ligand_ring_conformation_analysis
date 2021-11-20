#ifndef VECTOR_3D
#define VECTOR_3D

#include "point_3D.h"

struct Vector_3D : public Point_3D
{
        /* Constructors */
        Vector_3D(double _X = 0, double _Y = 0, double _Z = 0);
        Vector_3D(const Point_3D &A, const Point_3D &B = Point_3D(0, 0, 0));

        /* Arithmetic operations */
        double length() const;
        Vector_3D operator+(const Vector_3D &vec) const;
        Vector_3D operator-(const Vector_3D &vec) const;
        Vector_3D operator/(const double num) const;
        Vector_3D operator*(const double num) const;
        Vector_3D get_normal() const;
        static double dot(const Vector_3D &vec1, const Vector_3D &vec2);
        static Vector_3D cross(const Vector_3D &vec1, const Vector_3D &vec2);
};

#endif
