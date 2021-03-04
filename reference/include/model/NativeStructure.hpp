#pragma once
#include <utils/Types.hpp>
#include <istream>
#include <vector>

namespace CG {
    /* This class stores the native structure of a (part of a) chain. */
    class NativeStructure {
    public:
        int offset;

        /* Intra-chain native contacts. */
        struct Contact {
            std::pair<Index, Index> residues;
            Real distance;
        };
        std::vector<Contact> contacts;

        /* Tether and angle data:
         * bond(i) -> bond angle between i-1, i, i+1;
         * dihedral(i) -> dihedral angle between i-2, i-1, i, i+1;
         * tether(i) -> equilibrium distance between i, i+1. */
        RealList bond, dihedral, tether;
    };
}
