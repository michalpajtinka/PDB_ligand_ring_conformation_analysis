#ifndef SIX_ATOM_RING_H
#define SIX_ATOM_RING_H

#include "ring.h"
#include <string>

class Six_atom_ring: public Ring 
{
        public:
                Six_atom_ring() = delete;
                Six_atom_ring(std::string _structure);
                virtual ~Six_atom_ring() = 0;
        protected:
                /* functions for analyzing */
                virtual bool find_plane(double tolerance, int dist1 = 1, int dist2 = 3, int dist3 = 4);
                /* atom coordinates */
                Atom *C[6];
};

#endif
