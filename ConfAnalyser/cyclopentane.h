#ifndef CYCLOPENTANE_H
#define CYCLOPENTANE_H

#include "five_atom_ring.h"
#include <map>
#include <vector>

class Cyclopentane: public Five_atom_ring
{
        public:
                Cyclopentane() = delete;
                Cyclopentane(std::string _structure);
                virtual ~Cyclopentane();
                virtual bool analyse();
                virtual bool initialize(const std::vector<Atom*> &atoms);
        private:
                /* functions for analyzing */
                bool is_flat() const;
                bool is_envelope() const;
                bool is_half_chair() const;
		/* function verifying atom names for current ligand type */
		bool is_valid_atom_name(const int atom_number,
		                        const std::string &name) const;
                /* Tolerances */
                static constexpr double tolerance_in = 0.10;
                static constexpr double tolerance_out = 0.60;
                static constexpr double tolerance_tw_out = 0.54;
                static constexpr double angle_tw_boat = 10.5;
                static constexpr double angle_tolerance = 1;

};

#endif
