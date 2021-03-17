#pragma once
#include "forces/Force.hpp"
#include "data/PseudoAtoms.hpp"
#include "data/NativeStructure.hpp"
#include "utils/Units.hpp"

namespace cg::reference {
    using namespace cg::toolkit;

    class HarmonicTethers: public Force {
    private:
        PseudoAtoms const* pseudoAtoms;
        NativeStructure const* ns;

    public:
        Real default_dist0 = 3.8*angstrom;
        Real k1 = 100.0*eps/(angstrom*angstrom);
        Real k3 = 0.0;

        HarmonicTethers(PseudoAtoms const& pseudoAtoms, NativeStructure const& ns);

        void compute(Real &energy, Reals3 &forces) override;
    };
}