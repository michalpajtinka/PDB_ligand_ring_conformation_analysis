#include <algorithm>
#include <array>
#include <sstream>
#include <vector>
#include <string>
#include <cctype>
#include "oxane.h"
#include "plane_3D.h"
#include "angle.h"
#include "helper_functions.h"

#define ABOVE           0
#define UNDER           1

#define CHAIR           "CHAIR"
#define ENVELOPE        "ENVELOPE"
#define HALF_CHAIR      "HALF CHAIR"
#define BOAT            "BOAT"
#define SKEW            "SKEW"

using namespace std;


Oxane::Oxane(string _structure) : Six_atom_ring(_structure)
{
        conformations.insert({CHAIR, 3});
        conformations.insert({ENVELOPE, 4});
        conformations.insert({HALF_CHAIR, 5});
        conformations.insert({BOAT, 6});
        conformations.insert({SKEW, 7});

        for (auto x : outOfPlaneAtoms) {
                x.presence = false;
        }
}


Oxane::~Oxane()
{
        for (auto x : C) {
                if (x != nullptr) {
                        delete(x);
                }
        }
}


string Oxane::translate_conformation() const
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

        // TODO redo to switch case once conformtions uses ENUM values
        //      (switch cases labels must be compile time evaluable constant
        //      expressions)
        if (conformation == conformations[CHAIR]) {
                conf_name << "C";
        } else if (conformation == conformations[ENVELOPE]) {
                conf_name << "E";
        } else if (conformation == conformations[HALF_CHAIR]) {
                conf_name << "H";
        } else if (conformation == conformations[BOAT]) {
                conf_name << "B";
        } else if (conformation == conformations[SKEW]) {
                conf_name << "S";
        } else {
                for (auto conf : conformations) {
                        if (conf.second == conformation) {
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


bool Oxane::initialize(const vector<Atom*> &atoms)
{
        std::array<bool, 6> found {false};
        bool oxygen_found = false;

	for (auto x : atoms) {
	        if (ligand.empty()) {
                        ligand = x->get_residue_name();
                        if (atom_names.find(ligand) == atom_names.end()) {
                                cerr << "Ligand '" << ligand
                                        << "' not recognized!\n";
                                return false;
                        }
                }
                string tmp = strip(x->get_atom_name());
                for (int atom_index = 0; atom_index < 6; atom_index++) {
                        if (is_valid_atom_name(atom_index, tmp)) {
                                if (!filler(x, found[atom_index], C[atom_index])) {
                                        return false;
                                }
                                // identify oxygen atom
                                string element_name = strip(C[atom_index]->get_element_name());
                                if (element_name == "O") {
                                        if (oxygen_found) {
                                                cerr << "Oxygen atom found twice";
                                                return false;
                                        }
                                        oxygen_found = true;
                                        oxygen_position = atom_index;
                                }
                        }
                }
        }

        filled = all_of(found.begin(), found.end(), [](bool is_found){return is_found;});
        if (!filled) {
                cerr << "Not all atoms were found!" << endl;   
        }
        if (!oxygen_found) {
                cerr << "Unable to find oxygen atom position using element_name PDB field.\n";
        }

        return filled && oxygen_found;
}


bool Oxane::is_flat()
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


bool Oxane::is_chair()
{
        /* code below should be right in theory, but causes really big troubles
           to analysis...

        // atoms C1, C2, C4, and C5 has to lay in one plane
        // In fact, any four atom lying opposite each other have to lay in one plane
        // in CHAIR, but the symmetry of this conformation and numbering rules
        // stating that oxygen atom has to be nr. 6 force us to set begin to C1
        begin = oxygen_position % 6 + 1;  // start at C1
        Plane_3D plane(*(C[begin]), *(C[(begin+1)%6]), *(C[(begin+3)%6]));
        has_plane = abs(plane.distance_from(*(C[(begin+4)%6]))) < tolerance_in;
        if (!has_plane) {
                return false;
        }
        */

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
                outOfPlaneAtoms[0].position = right_dist > 0 ? ABOVE : UNDER;
                outOfPlaneAtoms[1].position = left_dist > 0 ? ABOVE : UNDER;
                outOfPlaneAtoms[0].atom_name = to_string(get_index_by_oxygen(2));
                outOfPlaneAtoms[1].atom_name = to_string(get_index_by_oxygen(5));
        }

        return isChair;
}


bool Oxane::is_half_chair()
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
                outOfPlaneAtoms[0].position = right_dist > 0 ? ABOVE : UNDER;
                outOfPlaneAtoms[1].position = left_dist > 0 ? ABOVE : UNDER;
                outOfPlaneAtoms[0].atom_name = to_string(get_index_by_oxygen(4));
                outOfPlaneAtoms[1].atom_name = to_string(get_index_by_oxygen(5));
        }

        return isHalfChair;
}


bool Oxane::is_boat()
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
                outOfPlaneAtoms[0].position = right_dist > 0 ? ABOVE : UNDER;
                outOfPlaneAtoms[1].position = left_dist > 0 ? ABOVE : UNDER;
                outOfPlaneAtoms[0].atom_name = to_string(get_index_by_oxygen(2));
                outOfPlaneAtoms[1].atom_name = to_string(get_index_by_oxygen(5));
        }

        /* sort outOfPlaneAtoms */
        sort(outOfPlaneAtoms, outOfPlaneAtoms+2, [](OutOfPlaneAtom x, OutOfPlaneAtom y){return x.atom_name < y.atom_name;});

        return isBoat;
}


bool Oxane::is_envelope()
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
                outOfPlaneAtoms[0].position = right_dist > 0 ? ABOVE : UNDER;
                outOfPlaneAtoms[1].position = left_dist > 0 ? ABOVE : UNDER;
                outOfPlaneAtoms[0].atom_name = to_string(get_index_by_oxygen(2));
                outOfPlaneAtoms[1].atom_name = to_string(get_index_by_oxygen(5));
        }

        return isEnv;
}


bool Oxane::is_skew()
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
                outOfPlaneAtoms[0].position = right_dist > 0 ? ABOVE : UNDER;
                outOfPlaneAtoms[1].position = left_dist > 0 ? ABOVE : UNDER;
                outOfPlaneAtoms[0].atom_name = to_string(get_index_by_oxygen(3));
                outOfPlaneAtoms[1].atom_name = to_string(get_index_by_oxygen(5));
        }

        return isSkew;
}


/* Convert to numbering relative to oxygen atom.
 * Oxygen atom is marked as 6, carbons will be marked 1 to 5 */
short unsigned Oxane::get_index_by_oxygen(short int delta_begin)  {
        int delta_oxygen = 6 - (oxygen_position);  // real O index vs. desired mark 6
        return ((begin + delta_begin   // atom number relative to `begin`
                + delta_oxygen         // shifted relative to oxygen position
                - 1)                   // -1 so that 6 stays 6 after modulo 6
                % 6 + 1);              // keep it in range 1 to 6
}


/* As oxane ring can be numbered so that oxygen atom is nr. 6 in both clockwise and
 * anticlockwise directions, every conformation name created based on
 * these numberings has at least one different descriptive name describing the
 * same conformation. For sake of clarity, use always the minimal possible
 * numbering (in case of ambiguity, prefer the clockwise solution).
 *
 * THIS FUNCTION DOES NOT CHECK WHETHER CONFORMATION DESCRIBED BY NAME IS
 * PHYSICALLY RELEVANT - E.G. EVEN THOUGH 4S2 AND 2S4 CONFORMATIONS ARE NOT LISTED
 * AMONG 6 POSSIBLE TYPES OF SKEW CONFORMATION IN THE LITERATURE, THE CONFORMATION
 * NAME WILL BE CORRECTLY GENERATED TO DESCRIBE WHAT IS IN THE PDB FILE */
void Oxane::fill_out_of_plane_atom_names(OutOfPlaneAtom &right_OOP_atom,
                                         short unsigned right_OOP_atom_delta_begin,
                                         OutOfPlaneAtom &left_OOP_atom,
                                         short unsigned left_OOP_atom_delta_begin) {
        // primary numbering
        array<int, 2> clockwise_positions {
                get_index_by_oxygen(right_OOP_atom_delta_begin),
                get_index_by_oxygen(left_OOP_atom_delta_begin)
        };

        // Create alternative numbering:
        // imagine this act as turning the oxane ring upside down in the only axis that
        // will keep the oxygen on the position with the number 6 - the one crossing
        // the central plane from the oxygen edge to the edge of the C3 atom.
        // Positions above/under will swap and regarding numbering, 4s and 5s will
        // change to 2s and 1s (and vice versa) while 3s and 6s stay unchanged.
        array<int, 2> anticlockwise_positions {
                6-(clockwise_positions[1])%6,
                6-(clockwise_positions[0])%6
        };

        // pick the one with lower avarage value
        if ((clockwise_positions[0] + clockwise_positions[1]) / 2.0
            <= (anticlockwise_positions[0] + anticlockwise_positions[1]) / 2.0) {
                right_OOP_atom.atom_name = to_string(clockwise_positions[0]);
                left_OOP_atom.atom_name = to_string(clockwise_positions[1]);
        } else {
                right_OOP_atom.atom_name = to_string(anticlockwise_positions[0]);
                left_OOP_atom.atom_name = to_string(anticlockwise_positions[1]);
        }
}


bool Oxane::analyse()
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


bool Oxane::is_valid_atom_name(const int atom_number,
                                        const std::string &name) const
{
        return (ligand.empty()) ? false :
                                any_of(atom_names[ligand][atom_number].begin(),
                                atom_names[ligand][atom_number].end(),
                                [&name](string &tmp){return tmp == name;});
}
