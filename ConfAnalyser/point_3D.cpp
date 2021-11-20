#include "point_3D.h"
#include <cmath>

using namespace std;

Point_3D::Point_3D(double X_arg, double Y_arg, double Z_arg)
{
        X = X_arg;
        Y = Y_arg;
        Z = Z_arg;
}


double Point_3D::distance_from(const Point_3D &A, const Point_3D &B)
{
        return hypot(hypot(A.X - B.X, A.Y - B.Y), A.Z - B.Z);
}
