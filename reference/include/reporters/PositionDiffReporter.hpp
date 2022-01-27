#pragma once
#include "reporters/Reporter.hpp"
#include "data/PseudoAtoms.hpp"
#include <fstream>
#include <map>

namespace cg::reference {
    class PositionDiffReporter: public Reporter {
    public:
        PositionDiffReporter(PseudoAtoms const &pseudoAtoms, string filePath);

        void report(int step) override;
    private:
        PseudoAtoms const &pseudoAtoms;
        map<unsigned, Reals3> referencePositions;
    };
}
