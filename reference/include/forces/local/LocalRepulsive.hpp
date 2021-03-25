#pragma once
#include "forces/Force.hpp"
#include "data/PseudoAtoms.hpp"
#include "data/NativeStructure.hpp"
#include "utils/Units.hpp"

namespace cg::reference {
    using namespace cg::toolkit;

    class LocalRepulsive: public Force {
    private:
        PseudoAtoms const& pseudoAtoms;
        NativeStructure const& ns;

    public:
        bool enabled = true;
        Real repulsive_cutoff = 5.0 * angstrom;
        Real force_cap = 1000.0;

        LocalRepulsive(PseudoAtoms const& pseudoAtoms, NativeStructure const& ns);

        void compute(Real &energy, Reals3 &forces) override;
    };
}
