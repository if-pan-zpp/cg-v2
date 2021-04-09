#include "forces/nonlocal/PauliExclusion.hpp"
using namespace cg::reference;
using namespace std;

PauliExclusion::PauliExclusion(PseudoAtoms const &pseudoAtoms,
                               Topology const &top,
                               Neighborhood const &verletList,
                               Real excludedRadius):
    top(top),
    pseudoAtoms(pseudoAtoms),
    verletList(verletList),
    excludedRadius(excludedRadius),
    sq_cutoff(excludedRadius * excludedRadius) {
    
}

void PauliExclusion::compute(Reals3 &forces) {
    energy = 0.;
    Reals3 const &pos = pseudoAtoms.pos;
    const Real sigma = excludedRadius / pow(2., 1./6.);

    for (auto const &[i, j] : verletList.pairs) {
        Real3 diff_vec = top.offset(pos.col(i), pos.col(j));
        Real sq_dist = diff_vec.squaredNorm();

        if (sq_dist < sq_cutoff) {
            Real dist = sqrt(sq_dist);
            Real rsi = sigma / dist;
            Real r6 = pow(rsi, 6); // TODO: check the assembly here
            energy += 4. * r6 * (r6 - 1.) + 1.;
            Real force = -24. * r6 * (1. - 2. * r6) / dist;

            if (force > force_cap) force = force_cap;
            if (force < -force_cap) force = -force_cap;

            force /= dist;
            forces.col(i) -= diff_vec * force;
            forces.col(j) += diff_vec * force;
        }
    }
}

void PauliExclusion::dumpResults(Results &results) {
    results.potEnergy += energy;
}
