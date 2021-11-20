#ifndef BENZENE_H
#define BENZENE_H

#include "six_atom_ring.h"
#include <map>
#include <vector>

class Benzene: public Six_atom_ring
{
        public:
                Benzene() = delete;
                Benzene(std::string _structure);
                virtual ~Benzene();
                virtual bool analyse();
                virtual bool initialize(const std::vector<Atom*> &atoms);
        private:
                /* functions for analyzing */
                bool is_flat() const;
                bool is_tw_boat() const;
		/* function verifying atom names for current ligand type */
		bool is_valid_atom_name(const int atom_number,
		                        const std::string &name) const;
                /* Tolerances */
                static constexpr double tolerance_flat_in = 0.1;
};

#endif
