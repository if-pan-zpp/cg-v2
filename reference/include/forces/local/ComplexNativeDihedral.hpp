#pragma once

#include "forces/Force.hpp"
#include "data/PseudoAtoms.hpp"
#include "data/NativeStructure.hpp"
#include "utils/Units.hpp"

namespace cg::reference {
    using namespace cg::toolkit;

    /* CPC14.pdf, 3.2.1, different for simple and complex
     * variants because of the hyperparameters. */
    class ComplexNativeDihedral : public Force {
    private:
        PseudoAtoms const *pseudoAtoms;
        NativeStructure const *ns;

    public:
        Real K1 = 0.33 * eps / (radian * radian);
        Real K3 = 0.33 * eps / (radian * radian);

        ComplexNativeDihedral(PseudoAtoms const &pseudoAtoms, NativeStructure const &ns);

        void compute(Real &energy, Reals3 &forces) override;
    };
}