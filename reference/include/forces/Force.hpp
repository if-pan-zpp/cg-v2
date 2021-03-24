#pragma once
#include "data/Primitives.hpp"

namespace cg::reference {
    class Force {
    public:
        /* Compute energy and force, add to buffers. */
        virtual void compute(Real& energy, Reals3& forces) = 0;
    };
}
