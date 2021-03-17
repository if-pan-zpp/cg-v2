#pragma once
#include "data/Primitives.hpp"

namespace cg::reference {
    class Integrator {
    public:
        virtual void step(int nsteps) = 0;
    };
}