#pragma once
#include "forces/Force.hpp"
#include "data/PseudoAtoms.hpp"
#include "data/NativeStructure.hpp"
#include "utils/Units.hpp"

namespace cg::reference {
    using namespace cg::toolkit;

    /* CPC14.pdf, 3.2 */
    class NativeBond: public Force {
    private:
        PseudoAtoms const* pseudoAtoms;
        NativeStructure const* ns;

    public:
        Real k = 30.0*eps/(radian*radian);

        NativeBond(PseudoAtoms const& pseudoAtoms, NativeStructure const& ns);

        void compute(Reals &energy, Reals3 &forces) override;
    };
}
