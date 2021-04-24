#pragma once
#include "forces/Force.hpp"
#include "data/PseudoAtoms.hpp"
#include "data/NativeStructure.hpp"
#include "utils/Units.hpp"

namespace cg::reference {
    using namespace cg::toolkit;

    class LocalRepulsive: public Force {
    private:
        PseudoAtoms const &pseudoAtoms;

    public:
        const Real repulsive_cutoff = 5.0 * angstrom;
        const Real sigma = repulsive_cutoff / pow(2., 1./6.);
        const Real force_cap = 200.0 * eps / angstrom;

        LocalRepulsive(PseudoAtoms const &pseudoAtoms);

        Real energy;

        void compute(Reals3 &forces) override;
        void dumpResults(Results &results) override;
    };
}
