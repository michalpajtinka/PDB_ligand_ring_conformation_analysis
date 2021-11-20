#include "plane_3D.h"
#include <cmath>


using namespace std;

Plane_3D::Plane_3D(double _a, double _b, double _c, double _d)
{
        a = _a;
        b = _b;
        c = _c;
        d = _d;
}


Plane_3D::Plane_3D(const Vector_3D &normal, double _d)
{
        init(normal, _d);
}


Plane_3D::Plane_3D(const Point_3D &A, const Point_3D &B, const Point_3D &C)
{
        Vector_3D _u = Vector_3D(B, A);
        Vector_3D _v = Vector_3D(B, C);
        Vector_3D normal = Vector_3D::cross(_u, _v);
        double _d = -(normal.X * A.X + normal.Y * A.Y + normal.Z * A.Z);
        init(normal, _d);
}


void Plane_3D::init(const Vector_3D &normal, double _d)
{
        a = normal.X;
        b = normal.Y;
        c = normal.Z;
        d = _d;
}



double Plane_3D::distance_from(const Point_3D &poi) const
{
        return (a * poi.X + b * poi.Y + c * poi.Z + d) / sqrt(a*a + b*b + c*c);
}


bool Plane_3D::is_on_plane(const Point_3D &poi, const double tolerance) const
{
        return abs(distance_from(poi)) <= tolerance;
}
