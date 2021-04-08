#pragma once
#include "forces/Force.hpp"
#include "data/PseudoAtoms.hpp"
#include "data/NativeStructure.hpp"
#include "utils/Units.hpp"

namespace cg::reference {
    using namespace cg::toolkit;

    class HarmonicTethers: public Force {
    public:
        Real H1 = 100.0*eps/(angstrom*angstrom);
        Real H2 = 0.0;
        Real force_cap = 1000.0 * f77unit;

        HarmonicTethers(PseudoAtoms const &pseudoAtoms, NativeStructure const &ns);

        void compute(Reals3 &forces) override;
        void dumpResults(Results &results) override;

    private:
        PseudoAtoms const &pseudoAtoms;
        NativeStructure const &ns;
        Real energy;
    };
}
