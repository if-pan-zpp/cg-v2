#pragma once
#include "data/Primitives.hpp"

namespace cg::reference {
    /* This is for native contacts and native bond/angle/tether data. */
    class NativeStructure {
    public:
        struct Contact {
            std::pair<unsigned, unsigned> residues;
            Real distance;
        };
        std::vector<Contact> contacts;
        Pairs ssbonds, nativeParts;
        Integers isNative;
        Reals bond, dihedral, tether;
        Pairs chains;
    };
}
