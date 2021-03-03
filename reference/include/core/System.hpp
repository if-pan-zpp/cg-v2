#pragma once
#include "utils/Types.hpp"

namespace CG {
    /* Particle system, i.e. residues, chains, SS bonds, native contacts,
     * bond and dihedral angles etc. */
    class System {
    public:
        System() = default;

        /* Residues data. For now everything common shall be stored here. */
        int nresidues;
        Real3List pos, vel;
        IndexPairList chain_ranges;
        IndexList chain_id;
    };
}