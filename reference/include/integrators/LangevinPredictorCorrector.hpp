#pragma once
#include <vector>
#include "integrators/Integrator.hpp"
#include "forces/Force.hpp"
#include "data/PseudoAtoms.hpp"
#include "utils/Units.hpp"

namespace cg::reference {
    using namespace cg::toolkit;

    class LangevinPredictorCorrector: public Integrator {
    public:
        LangevinPredictorCorrector(Real delta, PseudoAtoms &pseudoAtoms);
        void init(Reals3 &forces) override;
        void step(Reals3 &forces) override;

        static constexpr int K = 5; // order of predictor corrector method
        const Real gamma = 0.0; //noise turned off, original value is 2.0

    private:
        const Real delta;
        size_t n; // number of pseudoAtoms
        Reals3 *derivs[K + 1]; // derivatives of positions up to K'th order

        // derivs[0] and derivs[1] point to pseudoAtoms.pos and pseudoAtoms.vel.
        // All further derivatives are held here in highDerivatives.
        Reals3 highDerivatives[K - 1];

        const std::array<Real, K + 1> pred_corr_params =
        {3./16,  251./360,  1.,  11./18,  1./6,  1./60};
        const std::array<std::array<Real, K + 1>, K + 1> newton_symbol = 
            {{{{1, 1, 1, 1, 1,  1}},
            {{0, 1, 2, 3, 4,  5}},
            {{0, 0, 1, 3, 6, 10}},
            {{0, 0, 0, 1, 4, 10}},
            {{0, 0, 0, 0, 1,  5}},
            {{0, 0, 0, 0, 0,  1}}}};
    };
}
