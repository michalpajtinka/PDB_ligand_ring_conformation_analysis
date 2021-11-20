#include <cmath>
#include "vector_3D.h"

using namespace std ;


Vector_3D::Vector_3D(double _X, double _Y, double _Z) : Point_3D (_X, _Y, _Z) {}


Vector_3D::Vector_3D(const Point_3D &A, const Point_3D &B)
{
        X = A.X - B.X;
        Y = A.Y - B.Y;
        Z = A.Z - B.Z;
}


double Vector_3D::length() const
{
        return sqrt(pow(X, 2) + pow(Y, 2) + pow(Z, 2));
}


Vector_3D Vector_3D::operator+(const Vector_3D &vec) const
{
        return Vector_3D(X + vec.X, Y + vec.Y, Z + vec.Z);
}


Vector_3D Vector_3D::operator-(const Vector_3D &vec) const
{
        return Vector_3D(X - vec.X, Y - vec.Y, Z - vec.Z);
}


Vector_3D Vector_3D::operator/(const double num) const
{
        return Vector_3D(X / num, Y / num, Z / num);
}


Vector_3D Vector_3D::operator*(const double num) const
{
        return Vector_3D(X * num, Y * num, Z * num);
}


Vector_3D Vector_3D::get_normal() const
{
        float mag = length();

        return Vector_3D(X / mag, Y / mag, Z / mag);
}


double Vector_3D::dot(const Vector_3D &vec1, const Vector_3D &vec2)
{
        return vec1.X * vec2.X + vec1.Y * vec2.Y + vec1.Z * vec2.Z;
}


Vector_3D Vector_3D::cross(const Vector_3D &vec1, const Vector_3D &vec2)
{
        return Vector_3D(vec1.Y * vec2.Z - vec1.Z * vec2.Y,
                         vec1.Z * vec2.X - vec1.X * vec2.Z,
                         vec1.X * vec2.Y - vec1.Y * vec2.X);
}
