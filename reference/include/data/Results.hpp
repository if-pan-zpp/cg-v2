#pragma once
#include "data/Primitives.hpp"

namespace cg::reference {
    class Results {
    public:
        Real potEnergy;
        unsigned activeContacts;

        Results();
        void clear();
    };
};
