#pragma once
#include "data/Primitives.hpp"
#include "data/Results.hpp"

namespace cg::reference {
    class Force {
    public:
        /* Compute energy and force, add to buffers. */
        virtual void compute(Reals3& forces) = 0;
        virtual void dumpResults(Results &) = 0;
    };
}
