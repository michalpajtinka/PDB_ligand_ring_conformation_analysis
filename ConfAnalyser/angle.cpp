#define _USE_MATH_DEFINES

#include "angle.h"
#include <cmath>

static double rad_to_deg(double rad)
{
        return rad*180/M_PI;
}


double angle(const Vector_3D &vec1, const Vector_3D &vec2)
{
        double tmp = acos(Vector_3D::dot(vec1.get_normal(), vec2.get_normal()));
        return rad_to_deg(tmp);
}


double angle(const Point_3D &A, const Point_3D &B, const Point_3D &C)
{
        return angle(Vector_3D(B, A), Vector_3D(B, C));
}

double dihedral_angle(const Point_3D &A, const Point_3D &B,
                      const Point_3D &C, const Point_3D &D)
{
        Vector_3D _u = Vector_3D::cross(Vector_3D(B, A), Vector_3D(B, C));
        Vector_3D _v = Vector_3D::cross(Vector_3D(C, B), Vector_3D(C, D));
        Vector_3D normal = Vector_3D(B, C).get_normal();

        double tmp = atan2(Vector_3D::dot(Vector_3D::cross(_u, _v), normal),
                     Vector_3D::dot(_u, _v));
        return rad_to_deg(tmp);
}
