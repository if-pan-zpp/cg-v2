#pragma once
#include "forces/Force.hpp"
#include "data/PseudoAtoms.hpp"
#include "data/Topology.hpp"
#include "data/SharedData.hpp"
#include "utils/Units.hpp"
#include "Simulation.hpp"

namespace cg::reference {

    class LJWall: public Force {

        public:
            LJWall(PseudoAtoms const &pseudoAtoms,
                   Topology const &topology,
                   Simulation const &simulation,
                   LJWallSharedData &lj_shared_data,
                   WallSharedData &wall_shared_data);

            void compute(Reals3 &forces) override;
            void dumpResults(Results &results) override;

        private:
            PseudoAtoms const &pseudoAtoms;
            Topology const &topology;
            Simulation const &simulation;
            LJWallSharedData &lj_shared_data;
            WallSharedData &wall_shared_data;

            // TODO passing parameters
            Real wall_dist_raw = 5.;
            Real wall_dist = wall_dist_raw * angstrom;
            Real wall_pot_coeff = 4.;

            Real energy;

            enum ZWall {
                UPPER = true,
                LOWER = false
            };
    };
}