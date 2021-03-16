#pragma once
#include "data/Primitives.hpp"

namespace cg {
    class Force {
    public:
        /* Compute energy and force, add to buffers;
         * if nullptr, then don't compute that part.
         * Note: for now we can ignore energy, I guess */
        virtual void compute(Reals *energy, Reals3 *force) = 0;
    };
}