#include "forces/local/LocalRepulsive.hpp"
using namespace cg::reference;

LocalRepulsive::LocalRepulsive(PseudoAtoms const& _pseudoAtoms,
                                 NativeStructure const& _ns):
    pseudoAtoms(_pseudoAtoms),
    ns(_ns) {
    
}

void LocalRepulsive::compute(Real &energy, Reals3 &forces) {
    Reals3 forces_diff = Reals3::Zero(pseudoAtoms.n, 3);
    const Reals3 &positions = pseudoAtoms.pos;
    const Integers &chainId = ns.chainId;
    size_t residues = pseudoAtoms.n;
    Real energy_diff;

    for(size_t i = 0; i < residues - 2; i++) {
        if(chainId[i] == chainId[i+1] && chainId[i+1] == chainId[i+2]) {
            Real3 diff_vec = positions.row(i + 2) - positions.row(i);
            Real dist = diff_vec.norm();
            Real sq_dist = dist * dist;

            if(sq_dist < repulsive_cutoff * repulsive_cutoff) {
                Real rsi = repulsive_cutoff / ((pow(2., 1./6.)) * dist);
                Real r6 = pow(rsi, 6.);
                Real energy = 4 * r6 * (r6 - 1) + 1;
                Real force = 24 * r6 * (1 - 2 * r6) / dist;
                if (force > force_cap) force = force_cap;
                if (force < -force_cap) force = -force_cap;
                force /= -dist;
                forces.row(i) -= diff_vec * force;
                forces.row(i+2) += diff_vec * force;
            }
        }
    }

    forces += forces_diff;
    energy += energy_diff;
}