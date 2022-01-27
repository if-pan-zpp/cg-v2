#pragma once
#include "data/Primitives.hpp"

namespace cg::reference {
    class Integrator {
    public:
        virtual void init(Reals3 &forces) = 0;
        virtual void step(Reals3 &forces) = 0;
    };
}
