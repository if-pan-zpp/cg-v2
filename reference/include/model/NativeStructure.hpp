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

        /* Bond, dihedral angles. */
        RealList bond, dihedral;
    };
}
