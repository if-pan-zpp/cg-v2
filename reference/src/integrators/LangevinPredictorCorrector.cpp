#include "integrators/LangevinPredictorCorrector.hpp"
#include <random>
using namespace std;
using namespace cg::reference;

LangevinPredictorCorrector::LangevinPredictorCorrector(Real delta, PseudoAtoms &pseudoAtoms)
    : delta(delta) { 
    n = pseudoAtoms.n;
    for (size_t i = 0; i < K - 1; ++i) {
        highDerivatives[i] = Reals3::Zero(3, n);
    }
    
    derivs[0] = &pseudoAtoms.pos;
    derivs[1] = &pseudoAtoms.vel;
    for (size_t i = 0; i < K - 1; ++i) {
        derivs[i + 2] = highDerivatives + i;
    }
}

Real normalDistribution() {
    // TODO: add possibility to set seed
    // TODO: add rand from cg.f for comparison
    static mt19937 rng(0xC0FFEE);
    return normal_distribution<double>(0.0, 1.0)(rng);
}

void LangevinPredictorCorrector::init(Reals3 &forces) {
    *derivs[2] = forces * (0.5 * delta * delta);
}

void LangevinPredictorCorrector::step(Reals3 &forces) {
    Real temperature = 0.35; //TODO: get temperature from environment

    // Langevine dynamics
    const Real gamma2 = gamma / delta;
    const Real const2 = sqrt(2 * temperature * gamma * delta) * delta;

    Reals3 noise = Reals3::Zero(3, n);
    for (size_t i = 0; i < 3 * n; ++i) {
        noise(i) = normalDistribution();
    }

    *derivs[1] += const2 * noise;
    forces -= gamma2 * (*derivs[1]);
    
    // corrector:
    const Real deltsq = 0.5 * delta * delta;

    const Reals3 err = *derivs[2] - deltsq * forces;

    for (int der_nr = 0; der_nr <= K; der_nr++) {
        *derivs[der_nr] -= err * pred_corr_params[der_nr];
    }

    // predictor:
    for(int der_nr = 0; der_nr <= K; der_nr++) {
        for(int next_der = der_nr + 1; next_der <= K; next_der++) {
            *derivs[der_nr] += newton_symbol[der_nr][next_der] * *derivs[next_der];
        }
    }
}
