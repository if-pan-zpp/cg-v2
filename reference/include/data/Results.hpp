#pragma once
#include "data/Primitives.hpp"
#include "forces/nonlocal/QuasiAdiabatic.hpp"

namespace cg::reference {
    class Results {
    public:
        Real potEnergy;
        unsigned activeContacts;
        std::vector<QuasiAdiabatic::Contact> *qaContacts;

        Results();
        void clear();
    };
};
