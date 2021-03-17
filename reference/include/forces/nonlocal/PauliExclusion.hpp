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
    private:
        Topology const* top;
        Neighborhood const* excludedPairs;
        Real excludedRadius;

    public:
        PauliExclusion(Topology const& top,
            Real excludedRadius = 5.0*angstrom);

        void compute(Real &energy, Reals3 &forces) override;
    };
}