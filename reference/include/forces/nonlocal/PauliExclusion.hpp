#pragma once
#include "forces/Force.hpp"
#include "data/Topology.hpp"
#include "utils/Units.hpp"

namespace cg::reference {
    using namespace cg::toolkit;

    /* CPC14.pdf, 4.1 (repulsive part)
     * As per Łukasz's comment, we may compute these forces
     * for all relevant pairs, as opposed to only those not in
     * contact. */
    class PauliExclusion: public Force {
    public:
        PauliExclusion(PseudoAtoms const &pseudoAtoms,
                       Topology const &top,
                       Neighborhood const &verletList,
                       Real excludedRadius = 5.0*angstrom);

        void compute(Reals3 &forces) override;
        void dumpResults(Results &results) override;

    private:
        Topology const &top;
        PseudoAtoms const &pseudoAtoms;
        Neighborhood const &verletList;

        const Real excludedRadius;
        const Real sq_cutoff; // = excludedRadius ^ 2
        const Real force_cap = 200.0 * eps / angstrom;

        Real energy;
    };
}
