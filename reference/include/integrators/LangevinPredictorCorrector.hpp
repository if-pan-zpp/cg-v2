#pragma once
#include <vector>
#include "integrators/Integrator.hpp"
#include "forces/Force.hpp"
#include "data/PseudoAtoms.hpp"
#include "utils/Units.hpp"

namespace cg::reference {
    using namespace cg::toolkit;

    class LangevinPredictorCorrector: public Integrator {
    private:
        Real energy;
        Reals3 atomForces;
        std::vector<Force*> forceObjects;

    public:
        Real dt = 0.005*nanosecond;
        Real temperature; // What value?
        LangevinPredictorCorrector();

        void attachForce(Force *force);
        void step(int nsteps) override;
    };
}