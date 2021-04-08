#pragma once
#include "forces/Force.hpp"
#include "data/PseudoAtoms.hpp"
#include "data/NativeStructure.hpp"

namespace cg::reference {
    /* CPC14.pdf, 3.1 */
    class Chirality: public Force {
    private:
        PseudoAtoms const *pseudoAtoms;
        NativeStructure const *ns;

    public:
        Real echi = 1.0;
        Chirality(PseudoAtoms const &pseudoAtoms, NativeStructure const &ns);

        void compute(Real &energy, Reals3 &forces) override;
    };
}
