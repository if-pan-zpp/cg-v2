#pragma once
#include "integrators/Integrator.hpp"
#include "data/PseudoAtoms.hpp"
#include "utils/Units.hpp"

namespace cg::reference {
    using namespace cg::toolkit;

    class PredictorCorrector: public Integrator {
    public:
        PredictorCorrector(Real dt = 0.005*nanosecond);

        void step(int nsteps) override;
    };
}