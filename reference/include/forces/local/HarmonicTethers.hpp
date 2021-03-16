#pragma once
#include "forces/Force.hpp"
#include "data/System.hpp"
#include "data/NativeStructure.hpp"
#include "data/LocalStructure.hpp"
#include "utils/Units.hpp"

namespace cg {
    class HarmonicTethers: public Force {
    public:
        System *system;
        NativeStructure *ns;
        LocalStructure *ls;
        Real default_dist0 = 3.8*angstrom;
        Real k1 = 100.0*eps/(angstrom*angstrom);
        Real k3 = 0.0;

        void compute(Reals *energy, Reals3 *force) override;
    };
}