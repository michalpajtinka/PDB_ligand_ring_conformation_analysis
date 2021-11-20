#ifndef PLANE_3D_H
#define PLANE_3D_H

#include "vector_3D.h"


class Plane_3D
{
        public:
                Plane_3D(double _a, double _b, double _c, double _d);
                Plane_3D(const Vector_3D &normal, double _d);
                Plane_3D(const Point_3D &A,
                         const Point_3D &B,
                         const Point_3D &C);
                double distance_from(const Point_3D &poi) const;
                bool is_on_plane(const Point_3D &poi, double tolerance) const;
        private:
                void init(const Vector_3D &normal, const double _d);
                double a, b, c, d;
};


#endif
