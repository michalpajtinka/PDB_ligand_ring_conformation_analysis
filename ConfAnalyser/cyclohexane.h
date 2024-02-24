#ifndef CYCLOHEXANE_H
#define CYCLOHEXANE_H

#include "six_atom_ring.h"
#include <map>
#include <vector>

class Cyclohexane: public Six_atom_ring
{
        public:
                Cyclohexane() = delete;
                Cyclohexane(std::string _structure);
                virtual ~Cyclohexane();
                virtual bool analyse();
                virtual bool initialize(const std::vector<Atom*> &atoms);
        private:
                /* functions for analyzing */
                bool is_flat() const;
                bool is_half_chair() const;
                bool is_chair() const;
                bool is_boat() const; 
                bool is_tw_boat_right() const;
                bool is_tw_boat_left() const;
		/* function verifying atom names for current ligand type */
		bool is_valid_atom_name(const int atom_number,
		                        const std::string &name) const;
                /* Tolerances */
                static constexpr double tolerance_in = 0.1;
                static constexpr double tolerance_flat_in = 0.1;
                static constexpr double tolerance_out = 0.6;
                static constexpr double tolerance_tw_out = 0.4;
                static constexpr double angle_tw_boat = 17.1;
                static constexpr double angle_tolerance = 1;
};

#endif
