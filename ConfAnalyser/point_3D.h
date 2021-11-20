#ifndef POINT_3D_H
#define POINT_3D_H

struct Point_3D
{
        /* Constructor */
        Point_3D(double X_arg = 0, double Y_arg = 0, double Z_arg = 0);
        /* Methodes */
        static double distance_from(const Point_3D &A, const Point_3D &B);
        /* Coordinates */
        double X,Y,Z;
};

#endif
