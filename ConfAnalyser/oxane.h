#ifndef OXANE_H
#define OXANE_H

#include "six_atom_ring.h"
#include <map>
#include <vector>

struct OutOfPlaneAtom {
        std::string atom_name;
        int position;
        bool presence;
};

class Oxane : public Six_atom_ring
{
        public:
                Oxane() = delete;
                Oxane(std::string _structure);
                virtual ~Oxane();
                virtual bool analyse();
                virtual bool initialize(const std::vector<Atom*> &atoms);
                virtual std::string translate_conformation() const override;
        private:
                /* functions for analyzing */
                bool is_flat();
                bool is_half_chair();
                bool is_chair();
                bool is_boat(); 
                bool is_envelope();
                bool is_skew();

		/* function verifying atom names for current ligand type */
		bool is_valid_atom_name(const int atom_number,
		                        const std::string &name) const;

                /* Tolerances */
                static constexpr double tolerance_in = 0.1;
                static constexpr double tolerance_out = 0.3;
                
                /* Info about atoms lying out of plane */
                OutOfPlaneAtom outOfPlaneAtoms[2];
                signed short oxygen_position;
                short unsigned get_index_by_oxygen(short int delta_from_begin);
                void fill_out_of_plane_atom_names(
                                OutOfPlaneAtom &right_OOP_atom,
                                short unsigned right_OOP_atom_delta_begin,
                                OutOfPlaneAtom &left_OOP_atom,
                                short unsigned left_OOP_atom_delta_begin);
};

#endif
