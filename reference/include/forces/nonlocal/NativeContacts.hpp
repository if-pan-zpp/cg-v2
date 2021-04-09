#pragma once
#include "forces/Force.hpp"
#include "data/NativeStructure.hpp"
#include "data/Topology.hpp"
#include "utils/Units.hpp"
using std::vector;

namespace cg::reference {
    using namespace cg::toolkit;

    /* CPC14.pdf, 4.1
     * Note: here we may forego neighborhood lists */
    class NativeContacts: public Force {
    public:
        NativeContacts(PseudoAtoms const &pseudoAtoms,
                       NativeStructure const &ns,
                       Topology const &top,
                       Real cutoff);

        void compute(Reals3 &forces) override;
        void dumpResults(Results &results) override;

    private:
        PseudoAtoms const &pseudoAtoms;
        Topology const &top;

        const Real depth = 1. * eps;
        const Real force_cap = 1000. * eps/f77unit;
        const Real sq_cutoff;

        struct Contact {
            unsigned i, j;
            Real sigma;
        };
        vector<Contact> contacts;

        Real energy;
    };
}

