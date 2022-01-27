#pragma once
#include "reporters/Reporter.hpp"
#include "data/PseudoAtoms.hpp"
#include "data/Results.hpp"

namespace cg::reference {
    class StateReporter: public Reporter {
    public:
        StateReporter(PseudoAtoms const &pseudoAtoms, Results const &results, Real delta);

        void report(int step) override;
        double calculateAsphericity();

    private:
        Real delta;
        PseudoAtoms const &pseudoAtoms;
        Results const &results;
    };
}
