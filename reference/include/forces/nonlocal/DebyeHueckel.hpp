#pragma once
#include "forces/Force.hpp"
#include "data/PseudoAtoms.hpp"
#include "data/Topology.hpp"
#include "utils/Units.hpp"

namespace cg::reference {
    using namespace cg::toolkit;

    /* CPC14.pdf, 4.2.3 */
    class DebyeHueckel: public Force {
    private:
        PseudoAtoms const* pseudoAtoms;
        Topology const* top;
        Neighborhood const* neighborhood;

    public:
        bool includeAttractiveInteractions; /* !lrepcoul */
        Real screeningLength; /* screend */
        // Need to figure out the `coul` variable from CPC14.pdf

        /* Not sure what's the cutoff. */
        DebyeHueckel(PseudoAtoms const& pseudoAtoms, Topology const& top,
            Real cutoff = 0.0*angstrom);

        void compute(Real &energy, Reals3 &forces) override;
    };
}