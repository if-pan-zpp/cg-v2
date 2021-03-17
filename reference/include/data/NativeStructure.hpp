#pragma once
#include "data/Primitives.hpp"

namespace cg::reference {
    /* This is for native contacts and native bond/angle/tether data. */
    class NativeStructure {
    public:
        Pairs contacts, ssbonds, nativeParts;
        Integers isNative;
        Reals bond, dihedral, tether;
        Integers chainId;
        Pairs chains;
    };
}
