#pragma once
#include "math/Vector.hpp"

namespace CG {
    class System {
    public:
        /* Positions and velocities of the residues. */
        int nresidues;
        Real3List pos, vel;
        IndexList chain_id;
    };
}