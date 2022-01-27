#include "integrators/LangevinPredictorCorrector.hpp"
#include <random>
using namespace std;
using namespace cg::reference;

LangevinPredictorCorrector::LangevinPredictorCorrector(Real delta,
                                                       PseudoAtoms &pseudoAtoms,
                                                       RNG &rng)
    : delta(delta), rng(rng) {

    n = pseudoAtoms.n;
    masses = pseudoAtoms.mass;
    for (size_t i = 0; i < K - 1; ++i) {
        highDerivatives[i] = Reals3::Zero(3, n);
    }
    
    derivs[0] = &pseudoAtoms.pos;
    derivs[1] = &pseudoAtoms.vel;
    for (size_t i = 0; i < K - 1; ++i) {
        derivs[i + 2] = highDerivatives + i;
    }
}

void LangevinPredictorCorrector::init(Reals3 &forces) {
    for (size_t i = 0; i < n; ++i) {
        forces.col(i) /= masses(i);
    }

    *derivs[2] = forces * (0.5 * delta * delta);
}

Reals3 genGaussianNoise(size_t n, RNG &rng) {
    Eigen::Matrix<Real, 3, -1, Eigen::RowMajor> noise = Reals3::Zero(3, n);
    for (size_t i = 0; i < 3 * n; ++i) {
        noise(i) = rng.normal();
    }
    return noise;
}

void LangevinPredictorCorrector::step(Reals3 &forces) {
    Real temperature = 0.35 * toolkit::epsDivkB; //TODO: get temperature from environment

    // Langevine dynamics
    const Real gamma2 = gamma / delta;
    const Real const2 = sqrt(2 * temperature * gamma * delta) * delta;

    auto noise = genGaussianNoise(n, rng);
    for (size_t i = 0; i < n; ++i) {
        noise.col(i) /= masses[i];
    }
    *derivs[1] += const2 * noise;
    
    forces -= gamma2 * (*derivs[1]);
    for (size_t i = 0; i < n; ++i) {
        forces.col(i) /= masses[i];
    }

    const Real deltsq = 0.5 * delta * delta;
    const Reals3 err = *derivs[2] - deltsq * forces;

    for (int der_nr = 0; der_nr <= K; der_nr++) {
        *derivs[der_nr] -= err * pred_corr_params[der_nr];
    }

    for(int der_nr = 0; der_nr <= K; der_nr++) {
        for(int next_der = der_nr + 1; next_der <= K; next_der++) {
            *derivs[der_nr] += newton_symbol[der_nr][next_der] * *derivs[next_der];
        }
    }
}
