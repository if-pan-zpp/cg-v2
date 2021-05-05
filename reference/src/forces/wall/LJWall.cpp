#include "forces/wall/LJWall.hpp"
#include "data/Results.hpp"
using namespace cg::reference;

LJWall::LJWall(PseudoAtoms const &pseudoAtoms,
                   Topology const &topology,
                   Simulation const &simulation,
                   LJWallSharedData &lj_shared_data,
                   WallSharedData &wall_shared_data):
    pseudoAtoms(pseudoAtoms),
    topology(topology),
    simulation(simulation),
    lj_shared_data(lj_shared_data),
    wall_shared_data(wall_shared_data) {}


void LJWall::compute(Reals3 &forces) {
    size_t residues = pseudoAtoms.n;
    Reals3 const &positions = pseudoAtoms.pos;
    pair<Real3, Real3> bounds = topology.getBoundaries();
    vector<unsigned> &adiabCoeff = lj_shared_data.adiabCoeff;
    map<size_t, bool> &connected_to_zwall = lj_shared_data.connected_to_zwall;
    Reals3 &pull_ref_pos = lj_shared_data.pull_ref_pos;
    pair<Real3, Real3> &wall_forces = wall_shared_data.wall_forces;

    size_t &upper_n = wall_shared_data.upper_n;
    size_t &lower_n = wall_shared_data.lower_n;
    size_t &old_upper_n = wall_shared_data.old_upper_n;
    size_t &old_lower_n = wall_shared_data.old_lower_n;
    old_upper_n = old_lower_n = 0;

    energy = 0.;

    for(size_t i = 0; i+1 < residues; i++) {
        Real z_pos = positions.col(i)[2];
        Real z_bottom = bounds.first[2];
        Real z_top = bounds.second[2];
        bool stick_up, stick_down;

        if(z_top - z_pos <= wall_dist)  stick_up = true;
        if(z_pos - z_bottom <= wall_dist)  stick_down = true;

        if(!stick_up && !stick_down)    continue;

        Real dist;
        
        if(stick_down)  dist = z_pos - z_bottom;
        if(stick_up)    dist = z_pos - z_top;
        
        Real rsi = wall_dist_raw / dist;
        Real r6 = pow(rsi, 6.);

        if(adiabCoeff[i] == 0) { // & !fcc  ?
            if(simulation.allowWallConnect()) {
                adiabCoeff[i] = 2;
                if(stick_up && (upper_n < residues / 2)) {
                    connected_to_zwall[i] = UPPER;
                    upper_n ++; 
                }
                if(stick_down && (lower_n < residues / 2)) {
                    connected_to_zwall[i] = LOWER;
                    lower_n ++; 
                }
                pull_ref_pos.col(i) = positions.col(i);
                pull_ref_pos.col(i)[2] = 0.;
            }
        } else {
            if(stick_up)    old_upper_n++ ;
            if(stick_down)  old_lower_n++ ;
        }

        energy += wall_pot_coeff * (4. * r6 * (r6 - 1.) + 1.);
        Real force = wall_pot_coeff * 24. * r6 * (1. - 2. * r6) / dist;

        forces.col(i)[2] -= force;
        if(stick_down)  wall_forces.first[2] -= force;
        if(stick_up)    wall_forces.second[2] += force;
    }
}

void LJWall::dumpResults(Results &results) {
    results.potEnergy += energy;
}