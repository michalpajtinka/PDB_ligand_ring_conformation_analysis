#ifndef FIVE_ATOM_RING_H
#define FIVE_ATOM_RING_H

#include "ring.h"
#include <string>

class Five_atom_ring: public Ring 
{
        public:
                Five_atom_ring() = delete;
                Five_atom_ring(std::string _structure);
                virtual ~Five_atom_ring() = 0;
        protected:
                /* functions for analyzing */
                virtual bool find_plane(double tolerance, int dist1 = 1, int dist2 = 2, int dist3 = 3);
                /* atom coordinates */
                Atom *C[5];
};

#endif
