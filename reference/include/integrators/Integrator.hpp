#pragma once
#include "data/Primitives.hpp"

namespace cg::reference {
    class Integrator {
    public:
        virtual void step(Real delta, Reals3 &forces) = 0;
    };
}
