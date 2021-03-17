#pragma once
#include "forces/Force.hpp"
#include "data/PseudoAtoms.hpp"
#include "data/NativeStructure.hpp"
#include "utils/Units.hpp"

namespace cg::reference {
    using namespace cg::toolkit;

    /* CPC14.pdf, 3.2.1, different for simple and complex
     * variants because of the hyperparameters. */
    class SimpleNativeDihedral: public Force {
    private:
        PseudoAtoms const *pseudoAtoms;
        NativeStructure const *ns;

    public:
        Real k = 3.33*eps/(radian*radian);

        SimpleNativeDihedral(PseudoAtoms const& pseudoAtoms, NativeStructure const& ns);

        void compute(Reals *energy, Reals3 *force) override;
    };
}