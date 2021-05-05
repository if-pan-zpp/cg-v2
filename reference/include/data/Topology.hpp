#pragma once
#include "data/Primitives.hpp"
#include "data/PseudoAtoms.hpp"
#include "data/NativeStructure.hpp"
#include <vector>
#include <functional>
#include <memory>

namespace cg::reference {
    class Topology;

    class Neighborhood {
    public:
        struct Spec {
            int minLocalOffset;
            bool include4;
            Pairs* exclusions;
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
        Reals3 startPos;
        Real3 startCell;
        Real maxCorrectDist;
    };

    class Topology {
    private:
        friend class Neighborhood;
        std::vector<std::unique_ptr<Neighborhood>> neighborhoods;

        PseudoAtoms const* pseudoAtoms;
        NativeStructure const* ns;

        Real3 cell, cell_inv; /* 0 means no PBC along that axis. */

    public:
        Topology(PseudoAtoms const &pseudoAtoms, NativeStructure const &ns);

        /* Get the neighborhood with specificed specs. */
        Neighborhood const &createNeighborhood(Neighborhood::Spec const &spec);

        void updateCellSize(Real3 newCell);

        /* This is the vector from p to the closest image of q under the PBC
         * topology. */
        Real3 offset(Real3 const &p, Real3 const &q) const;

        /* return boundaries of simulation box (lowest point, highest point)
           zdown, xup etc. from cg.f */
        pair<Real3, Real3> getBoundaries() const; //TODO

        /* get shear value
        */
        Real getShear() const; // TODO

        /* Check whether neighborhoods are still correct. If not, recompute
         * them. */
        void update();
        void forceUpdate();
    };
}
