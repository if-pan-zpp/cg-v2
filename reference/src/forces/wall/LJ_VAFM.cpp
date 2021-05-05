#include "forces/wall/LJ_VAFM.hpp"
#include "data/Results.hpp"
using namespace cg::reference;

LJ_VAFM::LJ_VAFM(PseudoAtoms const &pseudoAtoms,
                   Topology const &topology,
                   LJWallSharedData &lj_shared_data,
                   WallSharedData &wall_shared_data):
    pseudoAtoms(pseudoAtoms),
    topology(topology),
    lj_shared_data(lj_shared_data),
    wall_shared_data(wall_shared_data) {}


void LJ_VAFM::compute(Reals3 &forces) {
    size_t residues = pseudoAtoms.n;
    Reals3 const &positions = pseudoAtoms.pos;
    pair<Real3, Real3> bounds = topology.getBoundaries();
    vector<unsigned> &adiabCoeff = lj_shared_data.adiabCoeff;
    map<size_t, bool> &connected_to_zwall = lj_shared_data.connected_to_zwall;
    Reals3 &pull_ref_pos = lj_shared_data.pull_ref_pos;
    pair<Real3, Real3> &wall_forces = wall_shared_data.wall_forces;
    Real shear = topology.getShear();

    size_t &old_upper_n = wall_shared_data.old_upper_n;
    size_t &old_lower_n = wall_shared_data.old_lower_n;

    Real &fresist = wall_shared_data.fresist;
    Real &fresistperp = wall_shared_data.fresistperp;

    energy = 0.;

    for (auto it = connected_to_zwall.begin(); it != connected_to_zwall.end(); it++) {
        size_t i = it->first;
        bool wall = it->second;
        Real zpull, shearpull;
        Real z_bottom = bounds.first[2];
        Real z_top = bounds.second[2];

        if(wall == LOWER) {
            zpull = z_bottom + pull_ref_pos[2];
            shearpull = shear;
            old_lower_n ++;
        }

        if(wall == UPPER) {
            zpull = z_top + pull_ref_pos[2];
            shearpull = - shear;
            old_upper_n ++;
        }

        Real3 diff_vec;
        diff_vec = positions.col(i);
        diff_vec[0] -= pull_ref_pos[0] + shearpull;
        diff_vec[1] -= pull_ref_pos[1];
        diff_vec[2] -= zpull;
        Real dist = diff_vec.norm();

        Real rsi = wall_dist_raw / dist;
        Real r6 = pow(rsi, 6.);
        Real pull_coeff = pulling_coeff * adiabCoeff[i] / max_adiab;
        energy = pull_coeff * 4. * r6 * (r6 - 1.);
        Real force = pull_coeff * 24. * r6 * (1. - 2. * r6) / dist;

        if(dist > wall_dist_raw * 1.5) {
            if(adiabCoeff[i] > 0) {
                adiabCoeff[i] --;
            } else {
                connected_to_zwall.erase(i);
            }
        } else {
            if(adiabCoeff[i] < max_adiab)
                adiabCoeff[i] ++;
        }

        if(dist > min_dist) {
            force /= - dist;
            if(force > force_cap) force = force_cap;
            else if(force < -force_cap) force = -force_cap;
            Real3 rep_force = force * diff_vec;
            forces += rep_force;

            if(wall == UPPER) {
                fresist -= rep_force[2];
                fresistperp -= rep_force[0];
            }

            if(wall == LOWER) {
                fresist += rep_force[2];
                fresistperp += rep_force[0];
            }
        }
    }
}


void LJ_VAFM::dumpResults(Results &results) {
    results.potEnergy += energy;
}