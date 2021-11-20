#include "ring.h"

using namespace std;

Ring::Ring(string _structure) : Molecule(_structure)
{
        conformations.insert({"FLAT", 2});
        has_plane = false;
        begin = 0;
}
