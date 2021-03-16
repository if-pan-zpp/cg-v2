#pragma once
#include "data/Primitives.hpp"
#include <unordered_map>

namespace cg {
    class System {
    public:
        Reals3 pos, vel;
        Real t;
        Reals mass, mass_inv;
        Names type;
        std::unordered_map<std::string, std::pair<int, int>> typeRanges;
    };
}