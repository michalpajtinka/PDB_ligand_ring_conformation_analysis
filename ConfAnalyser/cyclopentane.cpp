#include <algorithm>
#include "cyclopentane.h"
#include "plane_3D.h"
#include "angle.h"
#include "helper_functions.h"

#define ATOM_C1 0
#define ATOM_C2 1
#define ATOM_C3 2
#define ATOM_C4 3
#define ATOM_C5 4

using namespace std;


Cyclopentane::Cyclopentane(string _structure) : Five_atom_ring(_structure)
{
        conformations.insert({"ENVELOPE", 3});
        conformations.insert({"TWIST", 4});
}


Cyclopentane::~Cyclopentane()
{
        for (auto x : C) {
                if (x != nullptr) {
                        delete(x);
                }
        }
}


static bool filler(Atom *x, bool &found, Atom *&C)
{
        if (found) {
                cerr << strip(x->get_atom_name()) << " atom found twice!\n";
                return false;
        }

        found = true;

        C = new Atom(*x);

        return true;
}


bool Cyclopentane::initialize(const vector<Atom*> &atoms)
{
        bool found[5] = {false};
	for (auto x : atoms) {
	        if (ligand.empty()) {
                        ligand = x->get_residue_name();
                        if (atom_names.find(ligand) == atom_names.end()) {
                                cerr << "Ligand not recognized!" << endl;
                                return false;
                        }
                }
                string tmp = strip(x->get_atom_name());
                if (is_valid_atom_name(ATOM_C1, tmp)) {
                        if (!filler(x, found[0], C[0])) {
                                return false;
                        }
                } else if (is_valid_atom_name(ATOM_C2, tmp)) {
                        if (!filler(x, found[1], C[1])) {
                                return false;
                        } 
                } else if (is_valid_atom_name(ATOM_C3, tmp)) {
                        if (!filler(x, found[2], C[2])) {
                                return false;
                        } 
                } else if (is_valid_atom_name(ATOM_C4, tmp)) {
                        if (!filler(x, found[3], C[3])) {
                                return false;
                        } 
                } else if (is_valid_atom_name(ATOM_C5, tmp)) {
                        if (!filler(x, found[4], C[4])) {
                                return false;
                        } 
                }
        }

        filled = found[0] && found[1] && found[2] && found[3] && found[4];

        if (!filled) {
             cerr << "Not all atoms were found!" << endl;   
        }

        return filled;
}


bool Cyclopentane::is_flat() const
{
        /* flat conformation has all atoms in one plane */
        if (!has_plane) {
                return false;
        }

        Plane_3D plane(*(C[begin]), *(C[(begin+1)%5]), *(C[(begin+2)%5]));
        return plane.is_on_plane(*(C[(begin+3)%5]), tolerance_in) &&
               plane.is_on_plane(*(C[(begin+4)%5]), tolerance_in); 
}


bool Cyclopentane::is_envelope() const
{
        /* envelope conformation has all but one atom in one plane */
        if (!has_plane) {
                return false;
        }

        Plane_3D plane(*(C[begin]), *(C[(begin+1)%5]), *(C[(begin+2)%5]));

        return abs(plane.distance_from(*(C[(begin+4)%5]))) > tolerance_out;
}


bool Cyclopentane::is_twist() const
{
        /* twist conformation has no plane within the circle */
        if (has_plane) {
                return false;
        }

        Plane_3D left_plane(*(C[begin]),
                            *(C[(begin+1)%5]),
                            *(C[(begin+3)%5]));
        Plane_3D right_plane(*(C[begin]),
                             *(C[(begin+2)%5]),
                             *(C[(begin+3)%5]));
        double left_dist = left_plane.distance_from(*(C[(begin+4)%5]));
        double right_dist = right_plane.distance_from(*(C[(begin+4)%5]));
        double tw_angle = dihedral_angle(*(C[begin]), *(C[(begin+1)%5]),
					 *(C[(begin+2)%5]), *(C[(begin+3)%5]));
        return (abs(tw_angle) - angle_tw_boat < angle_tolerance) &&
               (abs(right_dist) > tolerance_tw_out) &&
               (abs(left_dist) > tolerance_tw_out) &&
               (right_dist * left_dist > 0);
}


bool Cyclopentane::analyse()
{
        if (!filled) {
                cerr << "Molecule has to be filled before analysis!" << endl;
                return false;
        }
        if (analysed) {
                cerr << "Attempt to analyze the same molecule twice!" << endl;
                return false;
        }

        /* finding the most accurate plane of 4 atoms within the ring */
        has_plane = find_plane(tolerance_in);
        if (is_flat()) {
                conformation = conformations["FLAT"];
        } else if (is_envelope()) {
                conformation = conformations["ENVELOPE"];
        } else if (is_twist()) {
		conformation = conformations["TWIST"];
        } else {
		conformation = conformations["UNDEFINIED"];
        }
                
        analysed = true;
        return true;
}


bool Cyclopentane::is_valid_atom_name(const int atom_number,
                                        const string &name) const
{
        return (ligand.empty()) ? false :
                                any_of(atom_names[ligand][atom_number].begin(),
                                atom_names[ligand][atom_number].end(),
                                [&name](string &tmp){return tmp == name;});
}
