#pragma once
#include "integrators/Integrator.hpp"
#include "data/System.hpp"

namespace cg {
    class PredictorCorrector: public Integrator {
    public:
        PredictorCorrector(Real dt);
        void step(int nsteps) override;
    };
}