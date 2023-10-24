#ifndef RING_H
#define RING_H

#include "molecule.h"
#include <string>

class Ring: public Molecule 
{
        public:
                Ring() = delete;
                Ring(std::string _structure);
        protected:
                /* Functions for analyzing */
                virtual bool find_plane(double tolerance, int dist1, int dist2, int dist3) = 0;

                /* Is the plane there? */
                bool has_plane;

                /* First atom of the plane */
                int begin;
};

#endif
