#include "five_atom_ring.h"
#include "plane_3D.h"
#include <cfloat>
#include <cmath>

using namespace std;

Five_atom_ring::Five_atom_ring(string _structure) : Ring(_structure)
{
        for (auto &x : C) {
                x = nullptr;
        }
}


Five_atom_ring::~Five_atom_ring() {}


bool Five_atom_ring::find_plane(double tolerance, int dist1, int dist2, int dist3)
{
        bool has_plane = false;
        double distance = DBL_MAX;
        for (int i = 0; i < 5; i++) {
                Plane_3D tmp(*(C[i%5]), *(C[(i+dist1)%5]), *(C[(i+dist2)%5]));
                double _distance = abs(tmp.distance_from(*(C[(i+dist3)%5])));
                if (_distance < distance) {
                        begin = i;
                        distance = _distance;
                        if (tmp.is_on_plane(*(C[(i+dist3)%5]), tolerance)) {
                                has_plane = true;
                        }
                }
        }

        return has_plane;
}
