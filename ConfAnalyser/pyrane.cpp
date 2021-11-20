#include <algorithm>
#include <sstream>
#include "pyrane.h"
#include "plane_3D.h"
#include "angle.h"
#include "helper_functions.h"

#define ATOM_C1 0
#define ATOM_C2 1
#define ATOM_C3 2
#define ATOM_C4 3
#define ATOM_C5 4
#define ATOM_O  5

#define ABOVE   0
#define UNDER   1

using namespace std;


Pyrane::Pyrane(string _structure) : Six_atom_ring(_structure)
{
        conformations.insert({"CHAIR", 3});
        conformations.insert({"ENVELOPE", 4});
        conformations.insert({"HALF CHAIR", 5});
        conformations.insert({"BOAT", 6});
        conformations.insert({"SKEW", 7});

        for (auto x : outOfPlaneAtoms) {
                x.presence = false;
        }
}


Pyrane::~Pyrane()
{
        for (auto x : C) {
                if (x != nullptr) {
                        delete(x);
                }
        }
}


string Pyrane::translate_conformation() const
{
        stringstream conf_name;

        bool first_symbol = true;

        for (auto x : outOfPlaneAtoms) {
                if (x.presence) {
                        if (x.position == ABOVE) {
                                if (!first_symbol) {
                                        conf_name << ",";
                                }
                                conf_name << x.atom_name;
                                first_symbol = false;
                        }
                }
        }

        for (auto conf : conformations) {
                if (conf.second == conformation) {
                        switch (conformation) {
                                case 3:
                                        conf_name << "C";
                                        break;
                                case 4:
                                        conf_name << "E";
                                        break;
                                case 5:
                                        conf_name << "H";
                                        break;
                                case 6:
                                        conf_name << "B";
                                        break;
                                case 7:
                                        conf_name << "S";
                                        break;
                                default:
                                        return conf.first;
                        }
                }
        }


        first_symbol = true;
        for (auto x : outOfPlaneAtoms) {
                if (x.presence) {
                        if (x.position == UNDER) {
                                if (!first_symbol) {
                                        conf_name << ",";
                                }
                                conf_name << x.atom_name;
                                first_symbol = false;
                        }
                }
        }

        return conf_name.str();
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


bool Pyrane::initialize(const vector<Atom*> &atoms)
{
        bool found[6] = {false};
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
                } else if (is_valid_atom_name(ATOM_O, tmp)) {
                        if (!filler(x, found[5], C[5])) {
                                return false;
                        } 
                }
        }

        filled = found[0] && found[1] && found[2] && found[3]
                                                && found[4] && found[5];
        if (!filled) {
             cerr << "Not all atoms were found!" << endl;   
        }

        return filled;
}


bool Pyrane::is_flat()
{
        has_plane = find_plane(tolerance_in);
        if (!has_plane) {
                return false;
        }

        bool isFlat = false;

        Plane_3D left_plane(*(C[begin]), *(C[(begin+1)%6]), *(C[(begin+4)%6]));
        Plane_3D right_plane(*(C[begin]), *(C[(begin+1)%6]), *(C[(begin+3)%6]));

        isFlat = left_plane.is_on_plane(*(C[(begin+2)%6]), tolerance_in) &&
               left_plane.is_on_plane(*(C[(begin+5)%6]), tolerance_in) &&
               right_plane.is_on_plane(*(C[(begin+2)%6]), tolerance_in) &&
               right_plane.is_on_plane(*(C[(begin+5)%6]), tolerance_in);

        return isFlat;
}


bool Pyrane::is_chair()
{
        has_plane = find_plane(tolerance_in);
        if (!has_plane) {
                return false;
        }

        bool isChair = false;

        Plane_3D left_plane(*(C[begin]), *(C[(begin+1)%6]), *(C[(begin+4)%6]));
        Plane_3D right_plane(*(C[begin]), *(C[(begin+1)%6]), *(C[(begin+3)%6]));
        double right_dist = (abs(right_plane.distance_from(*(C[(begin+2)%6]))) <
                             abs(left_plane.distance_from(*(C[(begin+2)%6])))) ?
                                right_plane.distance_from(*(C[(begin+2)%6])) :
                                left_plane.distance_from(*(C[(begin+2)%6]));
        double left_dist = (abs(right_plane.distance_from(*(C[(begin+5)%6]))) <
                             abs(left_plane.distance_from(*(C[(begin+5)%6])))) ?
                                right_plane.distance_from(*(C[(begin+5)%6])) :
                                left_plane.distance_from(*(C[(begin+5)%6]));

        isChair = (abs(right_dist) > tolerance_out &&
                abs(left_dist) > tolerance_out) &&
               (right_dist * left_dist < 0);

        if (isChair) {
                outOfPlaneAtoms[0].presence = true;
                outOfPlaneAtoms[1].presence = true;
                outOfPlaneAtoms[0].atom_name = ((begin+2)%6 == ATOM_O) ? "O" :
                                                     to_string((begin+2)%6 + 1);
                outOfPlaneAtoms[1].atom_name = ((begin+5)%6 == ATOM_O) ? "O" :
                                                     to_string((begin+5)%6 + 1); 
                outOfPlaneAtoms[0].position = right_dist > 0 ? ABOVE : UNDER;
                outOfPlaneAtoms[1].position = left_dist > 0 ? ABOVE : UNDER;
        }

        return isChair;
}


bool Pyrane::is_half_chair()
{
        has_plane = find_plane(tolerance_in, 1, 2, 3);
        if (!has_plane) {
                return false;
        }

        bool isHalfChair = false;

        Plane_3D left_plane(*(C[begin]), *(C[(begin+1)%6]), *(C[(begin+3)%6]));
        Plane_3D right_plane(*(C[begin]), *(C[(begin+1)%6]), *(C[(begin+2)%6]));
        double right_dist = (abs(right_plane.distance_from(*(C[(begin+4)%6]))) <
                             abs(left_plane.distance_from(*(C[(begin+4)%6])))) ?
                                right_plane.distance_from(*(C[(begin+4)%6])) :
                                left_plane.distance_from(*(C[(begin+4)%6]));
        double left_dist = (abs(right_plane.distance_from(*(C[(begin+5)%6]))) <
                             abs(left_plane.distance_from(*(C[(begin+5)%6])))) ?
                                right_plane.distance_from(*(C[(begin+5)%6])) :
                                left_plane.distance_from(*(C[(begin+5)%6]));

        isHalfChair = (abs(right_dist) > tolerance_out &&
                abs(left_dist) > tolerance_out) &&
               (right_dist * left_dist < 0);

        if (isHalfChair) {
                outOfPlaneAtoms[0].presence = true;
                outOfPlaneAtoms[1].presence = true;
                outOfPlaneAtoms[0].atom_name = ((begin+4)%6 == ATOM_O) ? "O" : to_string((begin+4)%6 + 1);
                outOfPlaneAtoms[1].atom_name = ((begin+5)%6 == ATOM_O) ? "O" : to_string((begin+5)%6 + 1); 
                outOfPlaneAtoms[0].position = right_dist > 0 ? ABOVE : UNDER;
                outOfPlaneAtoms[1].position = left_dist > 0 ? ABOVE : UNDER;
        }

        return isHalfChair;
}


bool Pyrane::is_boat()
{
        has_plane = find_plane(tolerance_in);
        if (!has_plane) {
                return false;
        }

        bool isBoat = false;

        Plane_3D left_plane(*(C[begin]), *(C[(begin+1)%6]), *(C[(begin+4)%6]));
        Plane_3D right_plane(*(C[begin]), *(C[(begin+1)%6]), *(C[(begin+3)%6]));
        double right_dist = (abs(right_plane.distance_from(*(C[(begin+2)%6]))) <
                             abs(left_plane.distance_from(*(C[(begin+2)%6])))) ?
                                right_plane.distance_from(*(C[(begin+2)%6])) :
                                left_plane.distance_from(*(C[(begin+2)%6]));
        double left_dist = (abs(right_plane.distance_from(*(C[(begin+5)%6]))) <
                             abs(left_plane.distance_from(*(C[(begin+5)%6])))) ?
                                right_plane.distance_from(*(C[(begin+5)%6])) :
                                left_plane.distance_from(*(C[(begin+5)%6]));

        isBoat = (abs(right_dist) > tolerance_out &&
                abs(left_dist) > tolerance_out) &&
               (right_dist * left_dist > 0);

        if (isBoat) {
                outOfPlaneAtoms[0].presence = true;
                outOfPlaneAtoms[1].presence = true;
                outOfPlaneAtoms[0].atom_name = ((begin+2)%6 == ATOM_O) ? "O" :
                                                     to_string((begin+2)%6 + 1);
                outOfPlaneAtoms[1].atom_name = ((begin+5)%6 == ATOM_O) ? "O" :
                                                     to_string((begin+5)%6 + 1);
                outOfPlaneAtoms[0].position = right_dist > 0 ? ABOVE : UNDER;
                outOfPlaneAtoms[1].position = left_dist > 0 ? ABOVE : UNDER;
        }

        /* sort outOfPlaneAtoms */
        sort(outOfPlaneAtoms, outOfPlaneAtoms+2, [](OutOfPlaneAtom x, OutOfPlaneAtom y){return x.atom_name < y.atom_name;});

        return isBoat;
}


bool Pyrane::is_envelope()
{
        has_plane = find_plane(tolerance_in);
        if (!has_plane) {
                return false;
        }

        bool isEnv = false;

        Plane_3D left_plane(*(C[begin]), *(C[(begin+1)%6]), *(C[(begin+4)%6]));
        Plane_3D right_plane(*(C[begin]), *(C[(begin+1)%6]), *(C[(begin+3)%6]));

        double right_dist = (abs(right_plane.distance_from(*(C[(begin+2)%6]))) <
                             abs(left_plane.distance_from(*(C[(begin+2)%6])))) ?
                                right_plane.distance_from(*(C[(begin+2)%6])) :
                                left_plane.distance_from(*(C[(begin+2)%6]));
        double left_dist = (abs(right_plane.distance_from(*(C[(begin+5)%6]))) <
                             abs(left_plane.distance_from(*(C[(begin+5)%6])))) ?
                                right_plane.distance_from(*(C[(begin+5)%6])) :
                                left_plane.distance_from(*(C[(begin+5)%6]));

        isEnv = ((left_plane.is_on_plane(*(C[(begin+2)%6]), tolerance_in) &&
                right_plane.is_on_plane(*(C[(begin+2)%6]), tolerance_in)) !=
                (left_plane.is_on_plane(*(C[(begin+5)%6]), tolerance_in) &&
                right_plane.is_on_plane(*(C[(begin+5)%6]), tolerance_in))) &&
                ((left_plane.is_on_plane(*(C[(begin+2)%6]), tolerance_in) ==
                right_plane.is_on_plane(*(C[(begin+2)%6]), tolerance_in)) &&
                (left_plane.is_on_plane(*(C[(begin+5)%6]), tolerance_in) ==
                right_plane.is_on_plane(*(C[(begin+5)%6]), tolerance_in)));

        if (isEnv) {
                outOfPlaneAtoms[0].presence = !left_plane.is_on_plane(*(C[(begin+2)%6]), tolerance_in);
                outOfPlaneAtoms[1].presence = !left_plane.is_on_plane(*(C[(begin+5)%6]), tolerance_in);
                outOfPlaneAtoms[0].atom_name = ((begin+2)%6 == ATOM_O) ? "O" :
                                                     to_string((begin+2)%6 + 1);
                outOfPlaneAtoms[1].atom_name = ((begin+5)%6 == ATOM_O) ? "O" :
                                                     to_string((begin+5)%6 + 1);
                outOfPlaneAtoms[0].position = right_dist > 0 ? ABOVE : UNDER;
                outOfPlaneAtoms[1].position = left_dist > 0 ? ABOVE : UNDER;
        }

        return isEnv;
}


bool Pyrane::is_skew()
{
        has_plane = find_plane(tolerance_in, 1, 2, 4);
        if (!has_plane) {
                return false;
        }

        bool isSkew = false;

        Plane_3D left_plane(*(C[begin]), *(C[(begin+1)%6]), *(C[(begin+4)%6]));
        Plane_3D right_plane(*(C[begin]), *(C[(begin+1)%6]), *(C[(begin+2)%6]));
        double right_dist = (abs(right_plane.distance_from(*(C[(begin+3)%6]))) <
                             abs(left_plane.distance_from(*(C[(begin+3)%6])))) ?
                                right_plane.distance_from(*(C[(begin+3)%6])) :
                                left_plane.distance_from(*(C[(begin+3)%6]));
        double left_dist = (abs(right_plane.distance_from(*(C[(begin+5)%6]))) <
                             abs(left_plane.distance_from(*(C[(begin+5)%6])))) ?
                                right_plane.distance_from(*(C[(begin+5)%6])) :
                                left_plane.distance_from(*(C[(begin+5)%6]));

        isSkew = (abs(right_dist) > tolerance_out &&
                abs(left_dist) > tolerance_out) &&
               (right_dist * left_dist < 0);

        if (isSkew) {
                outOfPlaneAtoms[0].presence = true;
                outOfPlaneAtoms[1].presence = true;
                outOfPlaneAtoms[0].atom_name = ((begin+3)%6 == ATOM_O) ? "O" : to_string((begin+3)%6 + 1);
                outOfPlaneAtoms[1].atom_name = ((begin+5)%6 == ATOM_O) ? "O" : to_string((begin+5)%6 + 1); 
                outOfPlaneAtoms[0].position = right_dist > 0 ? ABOVE : UNDER;
                outOfPlaneAtoms[1].position = left_dist > 0 ? ABOVE : UNDER;
        }

        return isSkew;
}


bool Pyrane::analyse()
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
        if (is_flat()) {
                conformation = conformations["FLAT"];
        } else if (is_chair()) {
                conformation = conformations["CHAIR"];
        } else if (is_half_chair()) {
                conformation = conformations["HALF CHAIR"];
        } else if (is_boat()) {
                conformation = conformations["BOAT"];
        } else if (is_envelope()) {
		conformation = conformations["ENVELOPE"];
        } else if (is_skew()) {
		conformation = conformations["SKEW"];
        } else {
		conformation = conformations["UNDEFINIED"];
        }
                
        analysed = true;
        return true;
}


bool Pyrane::is_valid_atom_name(const int atom_number,
                                        const std::string &name) const
{
        return (ligand.empty()) ? false :
                                any_of(atom_names[ligand][atom_number].begin(),
                                atom_names[ligand][atom_number].end(),
                                [&name](string &tmp){return tmp == name;});
}
