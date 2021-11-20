#ifndef ANGLE_H
#define ANGLE_H

#include "plane_3D.h"

double angle(const Vector_3D &vc1, const Vector_3D &vec2);
double angle(const Point_3D &A, const Point_3D &B, const Point_3D &C);
double dihedral_angle(const Point_3D &A, const Point_3D &B,
                      const Point_3D &C, const Point_3D &D);

#endif
