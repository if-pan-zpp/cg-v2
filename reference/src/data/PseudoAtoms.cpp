#include "data/PseudoAtoms.hpp"
using namespace cg::reference;
using cg::toolkit::RNG;

Real3 normalVector(RNG &rng) {
    Real3 res;

    for (size_t dim = 0; dim < 3; ++dim) {
    #ifndef LEGACY_MODE
        res(dim) = rng.normal();
    #else // cg.f:3807
        res(dim) = 2.0 * rng.uniform() - 1.0;
    #endif
    }

    return res.normalized();
}

void PseudoAtoms::initMovement(RNG &rng, Real temperature, Real delta) {
    #ifdef LEGACY_MODE
        // In cg.f:intvel3d residues' masses aren't considered, so we shadow them
        vector<Real> mass(n, 1.0);
    #endif

    for (size_t i = 0; i < n; ++i) {
        vel.col(i) = normalVector(rng);
    }
    
    // Shifting velocities so that total momentum is 0
    
    Real tot_mass = 0.0;
    Real3 tot_momentum = Real3::Zero();

    for (size_t i = 0; i < n; ++i) {
        tot_mass += mass[i];
        tot_momentum += mass[i] * vel;
    }
    for (size_t i = 0; i < n; ++i) {
        vel.col(i) -= tot_momentum * (1.0 / tot_mass);
    }

    // Scaling velocities to desired temperature

    Real average_kin_energy = 0.0;
    for (size_t i = 0; i < n; ++i) {
        average_kin_energy += mass[i] * vel.col(i).squaredNorm() / (2.0 * delta * delta);
    }
    average_kin_energy /= n;

    Real ratio = 1.5 * temperature / average_kin_energy;
    for (size_t i = 0; i < n; ++i) {
        vel.col(i) *= sqrt(ratio);
    }
}
