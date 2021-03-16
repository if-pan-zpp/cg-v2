#pragma once
#include "forces/Force.hpp"
#include "data/System.hpp"
#include "data/NativeStructure.hpp"
#include "data/LocalStructure.hpp"
#include "utils/Units.hpp"

namespace cg {
    class SimpleNativeDihedral: public Force {
    public:
        System *system;
        NativeStructure *ns;
        LocalStructure *ls;
        Real k = 3.33*eps/(radian*radian);

        void compute(Reals *energy, Reals3 *force) override;
    };
}