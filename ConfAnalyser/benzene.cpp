#include <algorithm>
#include "benzene.h"
#include "plane_3D.h"
#include "angle.h"
#include "helper_functions.h"

#define ATOM_C1 0
#define ATOM_C2 1
#define ATOM_C3 2
#define ATOM_C4 3
#define ATOM_C5 4
#define ATOM_C6 5

using namespace std;


Benzene::Benzene(string _structure) : Six_atom_ring(_structure) {}


Benzene::~Benzene()
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


bool Benzene::initialize(const vector<Atom*> &atoms)
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
                } else if (is_valid_atom_name(ATOM_C6, tmp)) {
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


bool Benzene::is_flat() const
{
        /* flat conforamtion has all atoms in one plane */
        if (!has_plane) {
                return false;
        }

        Plane_3D left_plane(*(C[begin]), *(C[(begin+1)%6]), *(C[(begin+4)%6]));
        Plane_3D right_plane(*(C[begin]), *(C[(begin+1)%6]), *(C[(begin+3)%6]));
        return left_plane.is_on_plane(*(C[(begin+2)%6]), tolerance_flat_in) &&
               left_plane.is_on_plane(*(C[(begin+5)%6]), tolerance_flat_in) &&
               right_plane.is_on_plane(*(C[(begin+2)%6]), tolerance_flat_in) &&
               right_plane.is_on_plane(*(C[(begin+5)%6]), tolerance_flat_in);
}


bool Benzene::analyse()
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
        has_plane = find_plane(tolerance_flat_in);
        if (is_flat()) {
                conformation = conformations["FLAT"];
        } else {
		conformation = conformations["UNDEFINIED"];
        }
                
        analysed = true;
        return true;
}


bool Benzene::is_valid_atom_name(const int atom_number,
                                        const std::string &name) const
{
        return (ligand.empty()) ? false :
                                any_of(atom_names[ligand][atom_number].begin(),
                                atom_names[ligand][atom_number].end(),
                                [&name](string &tmp){return tmp == name;});
}
