#pragma once
#include "forces/Force.hpp"
#include "data/NativeStructure.hpp"

namespace cg {
    class NativeContacts: public Force {
    public:
        NativeStructure *ns;

        void compute(Reals *energy, Reals3 *force) override;
    };
}

