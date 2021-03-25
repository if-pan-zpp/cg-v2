#pragma once
#include "reporters/Reporter.hpp"
#include "data/PseudoAtoms.hpp"

namespace cg::reference {
    class StateReporter: public Reporter {
    private:
        PseudoAtoms const& pseudoAtoms;

    public:
        StateReporter(PseudoAtoms const& pseudoAtoms);

        void report(int step) override;
    };
}
