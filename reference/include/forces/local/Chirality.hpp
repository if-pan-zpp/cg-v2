#pragma once
#include "forces/Force.hpp"
#include "data/System.hpp"
#include "data/NativeStructure.hpp"

namespace cg {
    class Chirality: public Force {
    private:
        System *system;
        NativeStructure *ns;
        Real echi;

    public:
        Chirality(System* system, NativeStructure *ns, Real echi);

        void compute(Reals *energy, Reals3 *force) override;
    };
}