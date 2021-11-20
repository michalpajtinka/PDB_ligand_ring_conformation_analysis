#include "six_atom_ring.h"
#include "plane_3D.h"
#include <cfloat>
#include <cmath>

using namespace std;

Six_atom_ring::Six_atom_ring(string _structure) : Ring(_structure)
{
        for (auto &x : C) {
                x = nullptr;
        }
}


Six_atom_ring::~Six_atom_ring() {}


// Old version
/*bool Six_atom_ring::find_plane(double tolerance)
{
        bool has_plane = false;
        double distance = DBL_MAX;
        for (int i = 0; i < 6; i++) {
                Plane_3D tmp(*(C[i%6]), *(C[(i+1)%6]), *(C[(i+3)%6]));
                double _distance = abs(tmp.distance_from(*(C[(i+4)%6])));
                if (_distance < distance) {
                        begin = i;
                        distance = _distance;
                        if (tmp.is_on_plane(*(C[(i+4)%6]), tolerance)) {
                                has_plane = true;
                        }
                }
        }

        return has_plane;
}*/


bool Six_atom_ring::find_plane(double tolerance, int dist1, int dist2, int dist3)
{
        bool has_plane = false;
        double distance = DBL_MAX;
        for (int i = 0; i < 6; i++) {
                Plane_3D tmp(*(C[i%6]), *(C[(i+dist1)%6]), *(C[(i+dist2)%6]));
                double _distance = abs(tmp.distance_from(*(C[(i+dist3)%6])));
                if (_distance < distance) {
                        begin = i;
                        distance = _distance;
                        if (tmp.is_on_plane(*(C[(i+dist3)%6]), tolerance)) {
                                has_plane = true;
                        }
                }
        }

        return has_plane;
}
