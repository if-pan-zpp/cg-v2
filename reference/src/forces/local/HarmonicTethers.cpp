#include "forces/local/HarmonicTethers.hpp"
using namespace cg::reference;

HarmonicTethers::HarmonicTethers(PseudoAtoms const& _pseudoAtoms,
                                 NativeStructure const& _ns):
    pseudoAtoms(_pseudoAtoms),
    ns(_ns) {
    
}

void HarmonicTethers::compute(Real &energy, Reals3 &forces) {
    size_t residues = pseudoAtoms.n;
    Reals3 forces_diff = Reals3::Zero(3, residues);
    Reals3 const& positions = pseudoAtoms.pos;
    Integers const& chainId = ns.chainId;
    Reals const& tether = ns.tether;
    Real energy_diff;

    for(size_t i = 0; i+1 < residues; i++) {
        if(chainId[i] == chainId[i+1]) {
            Real3 diff_vec = positions.col(i + 1) - positions.col(i);
            Real dist = diff_vec.norm();
            Real dist_change = dist - tether[i];
            Real sq_dist_change = dist_change * dist_change;
            Real energy_diff = (H1 + H2 * sq_dist_change) * sq_dist_change;
            Real force = (2 * H1 + 4 * H2 * sq_dist_change) * dist_change;

            if(force > force_cap) force = force_cap;
            else if(force < -force_cap) force = -force_cap;

            force /= -dist;
            forces_diff.col(i) -= diff_vec * force;
            forces_diff.col(i+1) += diff_vec * force;
        }
    }

    forces += forces_diff;
    energy += energy_diff;
}
