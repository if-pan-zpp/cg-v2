#pragma once

#include "forces/Force.hpp"
#include "data/PseudoAtoms.hpp"
#include "data/Topology.hpp"
#include "data/SharedData.hpp"
#include "data/Parameters.hpp"
#include "utils/Units.hpp"
using std::vector;

namespace cg::reference {
    class Results;
    
    class QuasiAdiabatic : public Force {
    public:

        QuasiAdiabatic(PseudoAtoms const &pseudoAtoms,
                       Topology const &top,
                       Neighborhood const &verletList,
                       toolkit::Parameters const &parameters,
                       SharedData &sharedData);

        void compute(Reals3 &forces) override;
        void dumpResults(Results &results) override;

        Real backbone_1_min = 0.75;
        Real backbone_2_min = 0.92;
        Real sidechain_max = 0.5;

        enum ContactType {
            NONE,
            BB,
            BS,
            SB,
            SS
        };
        bool disable(ContactType type, int i, int j);
        ContactType getContactType(int i, int j);

        struct Contact {
            unsigned i, j;
            ContactType type;
            unsigned adiabCoeff;
        };

    private:
        PseudoAtoms const &pseudoAtoms;
        Topology const &top;
        Neighborhood const &verletList;
        SharedData &sharedData;

        CoordNumber coordBounds[NUM_AMINO_ACIDS];

        enum AAType {
            GLY_OR_PRO,
            HYDROPHOBIC,
            POLAR,
            CHARGED_NEG,
            CHARGED_POS
        };
        AAType aaTypes[NUM_AMINO_ACIDS];

        uint8_t getSSBound(CoordNumber const &crdNum, AACode aaCode);

        Real3 hVector(int i);
        Real3 nVector(int i);

        vector<Contact> activeContacts;

        Real energy;
    };
}
