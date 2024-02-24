#include <algorithm>
#include "cyclohexane.h"
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


Cyclohexane::Cyclohexane(string _structure) : Six_atom_ring(_structure)
{
        conformations.insert({"CHAIR", 3});
        conformations.insert({"TW_BOAT_R", 4});
        conformations.insert({"TW_BOAT_L", 5});
        conformations.insert({"TW_BOAT_R", 4});
        conformations.insert({"TW_BOAT_L", 5});
        conformations.insert({"HALF_CHAIR", 6});
        conformations.insert({"BOAT", 7});
}


Cyclohexane::~Cyclohexane()
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


bool Cyclohexane::initialize(const vector<Atom*> &atoms)
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

/**
 * @brief Checks if the cyclohexane is in a flat conformation.
 * 
 * The analysis involves checking the relative positions of specific atoms in the cyclohexane ring.
 * - Flat conformation has all atoms in one plane
 * - The approach used in this function uses and analyzes two planes, right and left one, which should be almost identical.
 * - Left plane includes atoms 0, 1, 4.
 * - Right plane includes atoms 0, 1, 3.
 * - Atoms 2 and 5 should both lie on the left plane and on the right plane simultaneously
 * - The constant value <tolerance_flat_in> allows some deviation from that rule.
 * 
 * @return true 
 * @return false 
 */
bool Cyclohexane::is_flat() const
{
        if (!has_plane) {
                return false;
        }

        Plane_3D left_plane(*(C[begin]), *(C[(begin+1)%6]), *(C[(begin+4)%6]));
        Plane_3D right_plane(*(C[begin]), *(C[(begin+1)%6]), *(C[(begin+3)%6]));
        return left_plane.is_on_plane(*(C[(begin+2)%6]), tolerance_flat_in) &&
               left_plane.is_on_plane(*(C[(begin+5)%6]), tolerance_flat_in) &&  
               right_plane.is_on_plane(*(C[(begin+2)%6]), tolerance_flat_in) &&
               right_plane.is_on_plane(*(C[(begin+5)%6]), tolerance_flat_in);   // are atoms 2 and 5 on the right plane?
}

/**
 * @brief Checks if the cyclohexane is in a half-chair conformation.
 * 
 * The analysis involves checking the relative positions of specific atoms in the cyclohexane ring.
 * - Half-chair has all but one atom in one plane (atoms 0, 1, 3)
 * - Atom 2 or 5 (but not both) should not lie on the plane.
 * - Atom 4 should lie on the plane.
 * - Atom 2 or 5 should be far enough from the plane.
 * - The constants <tolerance_flat_in> and <tolerance_out> allow some deviations.
 * 
 * @return true 
 * @return false 
 */
bool Cyclohexane::is_half_chair() const
{
        if (!has_plane) {
                return false;
        }

        Plane_3D plane(*(C[begin]), *(C[(begin+1)%6]), *(C[(begin+3)%6]));
        double right_dist = plane.distance_from(*(C[(begin+2)%6]));
        double left_dist = plane.distance_from(*(C[(begin+5)%6]));
        return (plane.is_on_plane(*(C[(begin+2)%6]), tolerance_flat_in) !=
                plane.is_on_plane(*(C[(begin+5)%6]), tolerance_flat_in)) &&     
               plane.is_on_plane(*(C[(begin+4)%6]), tolerance_flat_in) &&       
               ((abs(right_dist) > tolerance_out) !=
                (abs(left_dist) > tolerance_out));                              
}

/**
 * @brief Checks if the cyclohexane is in a chair conformation.
 * 
 * The analysis involves checking the relative positions of specific atoms in the cyclohexane ring.
 * - Two atoms (2 and 5) of chair are on the opposite sides of plane (atoms 0, 1, 3)
 * - The atoms 2 and 5 should be far enough from the main plane. The constant <tolerance_out> allows some deviation.
 * 
 * @return true 
 * @return false 
 */
bool Cyclohexane::is_chair() const
{
        if (!has_plane) {
                return false;
        }
        Plane_3D plane(*(C[begin]), *(C[(begin+1)%6]), *(C[(begin+3)%6]));
        double right_dist = plane.distance_from(*(C[(begin+2)%6]));
        double left_dist = plane.distance_from(*(C[(begin+5)%6]));
        return (abs(right_dist) > tolerance_out &&
                abs(left_dist) > tolerance_out) &&
               (right_dist * left_dist < 0);            
}

/**
 * @brief Checks if the cyclohexane is in a boat conformation
 * 
 * - Two atoms (2 and 5) of boat are on the same side of a plane (atoms 0, 1, 3)
 * - The distance of atoms 2 and 5 should be far enough from the main plane.
 * - The constant <tolerance_out> allows some deviation.
 * 
 * @return true 
 * @return false 
 */
bool Cyclohexane::is_boat() const
{
        if (!has_plane) {
                return false;
        }
        Plane_3D plane(*(C[begin]), *(C[(begin+1)%6]), *(C[(begin+3)%6]));
        double right_dist = plane.distance_from(*(C[(begin+2)%6]));
        double left_dist = plane.distance_from(*(C[(begin+5)%6]));
        return (abs(right_dist) > tolerance_out &&
                abs(left_dist) > tolerance_out) &&      
               (right_dist * left_dist > 0);           
}

/**
 * @brief Checks if the cyclohexane is in a twisted boat conformation (boat with a right twist)
 * 
 * - Twisted boat has no plane within the ring
 * - The approach used in this function uses and analyzes two planes, right and left one
 * - Left plane includes atoms 0, 1, 4.
 * - Right plane includes atoms 0, 1, 3.
 * - The dihedral angle (1-3-4-0) should be in the allowed range (constant <angle_tolerance> allow some deviation)
 * - The atoms 2 and 5 (out of the plane) should be far enough from that plane. Constant <tolerance_tw_out> allows some deviation
 * - The atoms 2 and 5 should be on the same side (boat-like)
 * - The right plane should be further than the left one.
 * 
 * @return true 
 * @return false 
 */
bool Cyclohexane::is_tw_boat_right() const
{
        if (has_plane) {
                return false;
        }
        Plane_3D right_plane(*(C[begin]),
                             *(C[(begin+1)%6]),
                             *(C[(begin+3)%6]));
        Plane_3D left_plane(*(C[begin]),
                            *(C[(begin+1)%6]),
                            *(C[(begin+4)%6]));
        double right_dist = right_plane.distance_from(*(C[(begin+2)%6]));
        double left_dist = left_plane.distance_from(*(C[(begin+5)%6]));
        double tw_angle = dihedral_angle(*(C[(begin+1)%6]), *(C[(begin+3)%6]),
					 *(C[(begin+4)%6]), *(C[begin]));

        return ((abs(tw_angle) > angle_tw_boat - angle_tolerance) &&
                (abs(tw_angle) < angle_tw_boat + angle_tolerance)) &&   // is the twist angle in allowed range?
               ((abs(right_dist) > tolerance_tw_out) &&
                (abs(left_dist) > tolerance_tw_out)) &&                
                (right_dist * left_dist > 0) &&                  
                (right_dist > left_dist);                              
}

/**
 * @brief Checks if the cyclohexane is in a twisted boat conformation (boat with a left twist)
 * 
 * - Twisted boat has no plane within the ring
 * - The approach used in this function uses and analyzes two planes, right and left one
 * - Left plane includes atoms 0, 1, 4.
 * - Right plane includes atoms 0, 1, 3.
 * - The dihedral angle (1-3-4-0) should be in the allowed range (constant <angle_tolerance> allow some deviation)
 * - The atoms 2 and 5 (out of the plane) should be far enough from that plane. Constant <tolerance_tw_out> allows some deviation
 * - The atoms 2 and 5 should be on the same side (boat-like)
 * - The left plane should be further than the right one. (the only difference between this function and <is_tw_boat_right()>)
 * 
 * @return true 
 * @return false 
 */
bool Cyclohexane::is_tw_boat_left() const
{
        if (has_plane) {
                return false;
        }
        // Definition of plane made of atoms 0, 1 and 3
        Plane_3D right_plane(*(C[begin]),
                             *(C[(begin+1)%6]),
                             *(C[(begin+3)%6]));
        // Definition of plane made of atoms 0, 1 and 4
        Plane_3D left_plane(*(C[begin]),
                            *(C[(begin+1)%6]),
                            *(C[(begin+4)%6]));
        // Distance of atom 2 from right plane
        double right_dist = right_plane.distance_from(*(C[(begin+2)%6]));
        // Distance of atom 5 from left plane
        double left_dist = left_plane.distance_from(*(C[(begin+5)%6]));
        // Dihedral angle between atoms 1, 3, 4, 0
        double tw_angle = dihedral_angle(*(C[(begin+1)%6]), *(C[(begin+3)%6]),
					 *(C[(begin+4)%6]), *(C[begin]));

        return ((abs(tw_angle) > angle_tw_boat - angle_tolerance) &&
                (abs(tw_angle) < angle_tw_boat + angle_tolerance)) &&   
               ((abs(right_dist) > tolerance_tw_out) &&
                (abs(left_dist) > tolerance_tw_out)) &&                 
                (right_dist * left_dist > 0) &&                         
                (right_dist < left_dist);                              
}


bool Cyclohexane::analyse()
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
        } else if (is_half_chair()) {
                conformation = conformations["HALF_CHAIR"];
        } else if (is_boat()) {
                conformation = conformations["BOAT"];
        } else if (is_tw_boat_right()) {
		conformation = conformations["TW_BOAT_R"];
		conformation = conformations["TW_BOAT_R"];
        } else if (is_tw_boat_left()) {
		conformation = conformations["TW_BOAT_L"];
		conformation = conformations["TW_BOAT_L"];
        } else if (is_chair()) {
                conformation = conformations["CHAIR"];
        } else {
		conformation = conformations["UNDEFINIED"];
        }
                
        analysed = true;
        return true;
}


bool Cyclohexane::is_valid_atom_name(const int atom_number,
                                        const std::string &name) const
{
        return (ligand.empty()) ? false :
                                any_of(atom_names[ligand][atom_number].begin(),
                                atom_names[ligand][atom_number].end(),
                                [&name](string &tmp){return tmp == name;});
}
