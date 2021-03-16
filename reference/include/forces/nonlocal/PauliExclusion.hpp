#pragma once
#include "forces/Force.hpp"
#include "data/Contacts.hpp"
#include "utils/Units.hpp"

namespace cg {
    class PauliExclusion: public Force {
    public:
        Contacts *contacts;
        Real excludedRadius = 5.0*angstrom;

        void compute(Reals *energy, Reals3 *force) override;
    };
}