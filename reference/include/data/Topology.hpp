#pragma once
#include "data/Primitives.hpp"
#include "data/System.hpp"
#include "data/LocalStructure.hpp"
#include <vector>
#include <functional>

namespace cg {
    class Topology;

    class Neighborhood {
    public:
        struct Spec {
            int minLocalOffset;
            bool include4;
            Pairs* exclusions; /* must be lex-ordered */
            std::vector<std::pair<std::string, std::string>> pairTypes;
            Real cutoff, pad;
        };

        Neighborhood(Topology *owner, Spec spec);

    private:
        friend class Topology;
        Topology* top;
        Spec spec;

        void updateMaxCorrectDist();
        void update();
        void forceUpdate();

    public:
        Pairs pairs;
        Real3 startPos, startCell;
        Real maxCorrectDist;
    };

    class Topology {
    private:
        friend class Neighborhood;
        std::vector<Neighborhood> neighborhoods;

    public:
        System *system;
        LocalStructure *ls;

        Neighborhood const& createNeighborhood(Neighborhood::Spec const& spec);

        Real3 cell, cell_inv; /* 0 means no PBC along that axis. */
        void updateCellSize(Real3 newCell);

        /* This is the vector from p to the closest image of q under the PBC
         * topology. */
        Real3 offset(Real3 const& p, Real3 const& q) const;

        /* Check whether neighborhoods are still correct. If not, recompute
         * them. */
        void update();
        void forceUpdate();
    };
}