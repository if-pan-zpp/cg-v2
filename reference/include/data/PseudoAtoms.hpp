#pragma once
#include "data/Primitives.hpp"
#include <unordered_map>

namespace cg::reference {
    class PseudoAtoms {
    public:
        size_t n = 0; // number of pseudoatoms
        Reals3 pos, vel;
        Reals mass, mass_inv, charge;
        std::vector<std::string> type;
        std::unordered_map<std::string, std::pair<int, int>> typeRanges;
    };
}
