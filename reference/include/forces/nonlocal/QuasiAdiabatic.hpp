#pragma once

#include "forces/Force.hpp"
#include "data/PseudoAtoms.hpp"
#include "data/Topology.hpp"
#include "data/SharedData.hpp"
#include "utils/Units.hpp"

namespace cg::reference {
    using namespace cg::toolkit;
    using namespace std;

    class QuasiAdiabatic : public Force {
    public:
        enum ContactType {
            NONE,
            BB,
            BS,
            SB,
            SS
        };

        QuasiAdiabatic(PseudoAtoms const &pseudoAtoms,
                       Topology const &top,
                       SharedData &sharedData);

        void compute(Reals3 &forces) override;
        void dumpResults(Results &results) override;

        Real backbone_1_min = 0.75;
        Real backbone_2_min = 0.92;
        Real sidechain_max = 0.5;

        bool disable(ContactType type, int i, int j);
    private:
        PseudoAtoms const &pseudoAtoms;
        Topology const &top;
        SharedData &sharedData;


        Real3 hVector(int i);
        Real3 nVector(int i);
        ContactType getContactType(int i, int j);
    };
}
