#pragma once

#include "forces/Force.hpp"
#include "data/Topology.hpp"
#include "data/SharedData.hpp"
#include "utils/Units.hpp"

/*
  Technically, this isn't a force, but it fits into the interface.
  TODO: create separate interface, called 'Computation' or sth like that
 */

namespace cg::reference {
    class NeighborsCounter : public Force {
    public:
        NeighborsCounter(PseudoAtoms const &pseudoAtoms,
                         NativeStructure const &ns,
                         Topology const &top,
                         Neighborhood const &verletList,
                         SharedData &sharedData);

        void compute(Reals3 &) override;
        void dumpResults(Results &) override;

        const Real cutoff = 7.5 * toolkit::angstrom;
        
    private:
        PseudoAtoms const &pseudoAtoms;
        Topology const &top;
        Neighborhood const &verletList;
        SharedData &sharedData;

        vector<unsigned> neiCount;
        std::vector<pair<int, int>> nativeContacts;

        void checkPair(int i, int j);
        const Real cutoffSquared = cutoff * cutoff;
    };
}
