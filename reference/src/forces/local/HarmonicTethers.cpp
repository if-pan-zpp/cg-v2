#include "forces/local/HarmonicTethers.hpp"
using namespace cg::reference;

HarmonicTethers::HarmonicTethers(PseudoAtoms const& _pseudoAtoms,
                                 NativeStructure const& _ns):
    pseudoAtoms(_pseudoAtoms),
    ns(_ns) {
    
}

void HarmonicTethers::compute(Real &energy, Reals3 &forces) {
    Reals3 forces_diff = Reals3::Zero(pseudoAtoms.n, 3);
    const Reals3 &positions = pseudoAtoms.pos;
    const Integers &chainId = ns.chainId;
    size_t residues = pseudoAtoms.n;
    Reals tether = ns.tether;
    Real energy_diff;

    for(size_t i = 0; i < residues - 1; i++) {
        if(chainId[i] == chainId[i+1]) {
            Real3 diff_vec = positions.row(i + 1) - positions.row(i);
            Real dist = diff_vec.norm();
            Real dist_change = dist - tether[i];
            Real sq_dist_change = dist_change * dist_change;
            Real energy_diff = (H1 + H2 * sq_dist_change) * sq_dist_change;
            Real force = (2 * H1 + 4 * H2 * sq_dist_change) * dist_change;

            if(force > force_cap) force = force_cap;
            else if(force < -force_cap) force = -force_cap;

            force /= -dist;
            forces_diff.row(i) -= diff_vec * force;
            forces_diff.row(i+1) += diff_vec * force;
        }
    }

    forces += forces_diff;
    energy += energy_diff;
}
