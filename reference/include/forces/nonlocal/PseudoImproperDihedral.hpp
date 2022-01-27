#pragma once
#include "forces/Force.hpp"
#include "data/PseudoAtoms.hpp"
#include "data/NativeStructure.hpp"
#include "utils/Units.hpp"
#include "data/Topology.hpp"

namespace cg::reference {
    using namespace cg::toolkit;

    class PseudoImproperDihedral: public Force {
    private:
        PseudoAtoms const &pseudoAtoms;
        NativeStructure const &ns;
        Neighborhood const &verlet_list;

    public:
        bool enabled = false;
        Real cutoff = 18.0 * angstrom;
        bool use_mj_matrix = false;
        //PairMatrix mjMatrix; // TODO make the three commented things work
        //PairMatrix ssMinimalDistances;
        bool pid_cos = false;
        Real alpha_bb_pos = 6.4; // TODO add units
        Real alpha_bb_neg = 6.0;
        Real alpha_ss = 1.2;
        Real psi0_bb_pos = 1.05;
        Real psi0_bb_neg = -1.44;
        Real psi0_ss = -0.23;
        Real rmin_pos = 5.6 * angstrom;
        Real rmin_neg = 6.2 * angstrom;
        Real contact_mltp = 1.3;
        bool sink_pot = false;
        Real eps_bb = 0.2;
        bool pid_barrier = false;
        //electrostatics
        //Polarity polarity;
        bool pid_electrostatics = false;
        Real elektr_screen = 10.0 * angstrom;
        Real coul = 85.0 / (angstrom*angstrom);
        bool ele_perm_const = false;

        Real min_lambda = 0.00005;
        Real min_norm = 0.01;

        PseudoImproperDihedral(PseudoAtoms const &pseudoAtoms, NativeStructure const &ns, 
                                Neighborhood const &verlet_list);

        void compute(Reals3 &forces) override;
    };
}
