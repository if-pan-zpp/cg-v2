#include "forces/nonlocal/NativeContacts.hpp"
using namespace cg::reference;

NativeContacts::NativeContacts(PseudoAtoms const &pseudoAtoms,
                               NativeStructure const &ns,
                               Topology const &top,
                               Real cutoff):
    pseudoAtoms(pseudoAtoms),
    top(top),
    sq_cutoff(cutoff * cutoff) {

    Real scaling_constant = pow(0.5, 1./6.);
    for (NativeStructure::Contact const &contact : ns.contacts) {
        contacts.push_back(Contact {
                .i = contact.residues.first,
                .j = contact.residues.second,
                .sigma = contact.distance * scaling_constant
            });
    }
}

void NativeContacts::compute(Reals3 &forces) {
    Reals3 const &pos = pseudoAtoms.pos;
    energy = 0.;
    activeContacts = 0;

    for (Contact const &contact : contacts) {
        Reals3 diff_vec = top.offset(pos.col(contact.i), pos.col(contact.j));

        Real sq_dist = diff_vec.squaredNorm();
        if (sq_dist > sq_cutoff) continue;

        activeContacts++;

        Real dist = sqrt(sq_dist);
        Real sigma_by_dist_6 = pow(contact.sigma / dist, 6.); //TODO: check assembly

        energy += 4. * depth * sigma_by_dist_6 * (sigma_by_dist_6 - 1.);
        Real force = 24. * depth * sigma_by_dist_6 * (1. - 2. * sigma_by_dist_6) / dist;

        if (force > force_cap) force = force_cap;
        if (force < -force_cap) force = -force_cap;

        force /= dist;
        forces.col(contact.i) += diff_vec * force;
        forces.col(contact.j) -= diff_vec * force;
    }
}

void NativeContacts::dumpResults(Results &results) {
    results.potEnergy += energy;
    results.activeContacts += activeContacts;
}
