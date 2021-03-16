#pragma once
#include "data/Primitives.hpp"

namespace cg {
    class Integrator {
    public:
        virtual void step(int nsteps) = 0;
    };
}