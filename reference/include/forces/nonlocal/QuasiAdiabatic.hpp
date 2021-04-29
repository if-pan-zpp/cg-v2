#pragma once

#include "forces/Force.hpp"
#include "data/PseudoAtoms.hpp"
#include "data/Topology.hpp"
#include "data/SharedData.hpp"
#include "data/Parameters.hpp"
#include "utils/Units.hpp"
using std::vector;
using cg::toolkit::angstrom;

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
            BB,
            BS,
            SB,
            SS,
            NONE
        };
        static const size_t NUM_CONTACT_TYPES = ContactType::NONE;

        bool isDisabled(ContactType type, int i, int j) const;
        ContactType getContactType(int i, int j) const;

        Real getCutoff(AACode i, AACode j, ContactType type) const;

        struct Contact {
            unsigned i, j;
            ContactType type;
            unsigned adiabCoeff;
        };

        const Real cutoffCoeff = 1.0 * pow(2., 1./6.);

        const Real lj_r_min[NUM_CONTACT_TYPES] = {
            [BB] = 5.0 * angstrom,
            [BS] = 6.6 * angstrom,
            [SB] = 6.6 * angstrom,
            [SS] = 7.5 * angstrom
        };
        Real ss_lj_r_min[NUM_AMINO_ACIDS][NUM_AMINO_ACIDS];

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

        uint8_t getSSBound(CoordNumber const &crdNum, AACode aaCode) const;
        Real getMaxCutoff(AACode i, AACode j) const;

        Real3 hVector(int i) const;
        Real3 nVector(int i) const;

        vector<Contact> activeContacts;

        Real energy;

        Real totalCutoff, totalCutoffSquared;
    };
}
