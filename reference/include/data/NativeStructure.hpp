#pragma once
#include "data/Primitives.hpp"

namespace cg {
    class NativeStructure {
    public:
        Pairs contacts, ssbonds, nativeParts;
        Integers isNative;
        Reals bond, dihedral, tether;
    };
}
