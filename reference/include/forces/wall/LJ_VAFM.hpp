#pragma once
#include "forces/Force.hpp"
#include "data/PseudoAtoms.hpp"
#include "data/Topology.hpp"
#include "data/SharedData.hpp"
#include "utils/Units.hpp"

namespace cg::reference {

    class LJ_VAFM: public Force {

        public:
            LJ_VAFM(PseudoAtoms const &pseudoAtoms,
                   Topology const &topology,
                   LJWallSharedData &lj_shared_data,
                   WallSharedData &wall_shared_data);

            void compute(Reals3 &forces) override;
            void dumpResults(Results &results) override;

        private:
            PseudoAtoms const &pseudoAtoms;
            Topology const &topology;
            LJWallSharedData &lj_shared_data;
            WallSharedData &wall_shared_data;

            // TODO passing parameters
            unsigned max_adiab = 2000;
            Real force_cap = 1000.;
            Real wall_dist_raw = 5.;
            Real pulling_coeff = 1.;
            Real min_dist = 1e-10;
            Real energy;

            enum ZWall {
                UPPER = true,
                LOWER = false
            };
    };
}