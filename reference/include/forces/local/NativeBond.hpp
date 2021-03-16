#pragma once
#include "forces/Force.hpp"
#include "data/System.hpp"
#include "data/NativeStructure.hpp"
#include "utils/Units.hpp"

namespace cg {
    class NativeBond: public Force {
    public:
        System *system;
        NativeStructure *ns;
        Real k = 30.0*eps/(radian*radian);

        void compute(Reals *energy, Reals3 *force) override;
    };
}
