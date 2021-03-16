#pragma once
#include "forces/Force.hpp"
#include "data/System.hpp"
#include "data/NativeStructure.hpp"
#include "data/LocalStructure.hpp"
#include "utils/Units.hpp"

namespace cg {
    class ComplexNativeDihedral: public Force {
    public:
        System *system;
        NativeStructure *ns;
        LocalStructure *ls;
        Real K1 = 0.33*eps/(radian*radian);
        Real K3 = 0.33*eps/(radian*radian);

        void compute(Reals *energy, Reals3 *force) override;
    };
}