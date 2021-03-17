#pragma once
#include "data/Primitives.hpp"
#include <unordered_map>

namespace cg::reference {
    class PseudoAtoms {
    public:
        Reals3 pos, vel;
        Real time;
        Reals mass, mass_inv, charge;
        std::vector<std::string> type;
        std::unordered_map<std::string, std::pair<int, int>> typeRanges;
    };
}