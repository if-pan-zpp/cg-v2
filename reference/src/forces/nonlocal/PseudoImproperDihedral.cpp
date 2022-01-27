#include "forces/nonlocal/PseudoImproperDihedral.hpp"
using namespace cg::reference;

PseudoImproperDihedral::PseudoImproperDihedral(PseudoAtoms const &_pseudoAtoms,
                                 NativeStructure const &_ns, Neighborhood const &_verlet_list):
    pseudoAtoms(pseudoAtoms),
    ns(ns),
    verlet_list(verlet_list) {
}

void PseudoImproperDihedral::compute(Reals3 &forces) {
    size_t residues = pseudoAtoms.n;
    Reals3 const &positions = pseudoAtoms.pos;
    std::vector<std::string> const &type = pseudoAtoms.type;
    Pairs const &list = verlet_list.pairs;

    for (auto it = list.begin(); it != list.end(); it++) {
        int i = it -> first;
        int j = it -> second;

        bool norm_too_small = false;
        Real eps_mj;
        //TODO use MJ matrix if needed
        eps_mj = 1.;

        // TODO implement PBC distance
        Real3 diff_vec = positions.col(j) - positions.col(i);
        Real dist = diff_vec.norm();
        Real sq_dist = dist * dist;

        // TODO update residues' neighbour count

        if(sq_dist > cutoff * cutoff)   continue;

        if(type[i] == "PRO" || type[j] == "PRO")  continue;

        Real ss_lambda[2];
        Real bb_lambda[2];
        Real alpha1[2];
        Real alpha2[2];
        Real alpha3[2]; // a(_,4) in cg.f
        Real r_min;
        Real3 f_var[2][4]; // TODO find out what it represents and name it properly
        Real dvdp[2]; // TODO find descriptive name for it

        for(size_t nr = 0; nr <=1; nr++) {
            int i1 = i;
            int i2 = j;
            if(nr == 1) std::swap(i1,i2);
            // TODO check if i1-1, i1+1 etc. may be out of bounds
            Real3 v1 = positions.col(i1) - positions.col(i1+1);
            Real3 v2 = positions.col(i1-1) - positions.col(i1+1);
            Real3 v3 = positions.col(i1-1) - positions.col(i2);
            Real3 v4 = v1.cross(v2);
            Real3 v5 = v2.cross(v3);
            Real v4_norm_sq = v4.dot(v4);
            Real v5_norm_sq = v5.dot(v5);

            if(v4_norm_sq < min_norm || v5_norm_sq < min_norm) {
                norm_too_small = true;
                continue;
            }

            Real cospsi = v5.dot(v4) / sqrt(v5_norm_sq * v4_norm_sq);
            Real psi = acos(cospsi);
            if(v1.dot(v5) < 0)  psi *= -1.;

            alpha1[nr] = alpha_ss * (psi - psi0_ss);
            if(alpha1[nr] < M_PI && alpha1[nr] > -M_PI) {
                if(pid_cos) {
                    ss_lambda[nr] = 0.5 * (cos(alpha1[nr]) + 1);
                } else {
                    Real alpha_sq = (alpha1[nr] / M_PI) * (alpha1[nr] / M_PI);
                    ss_lambda[nr] = 1. - alpha_sq 
                                    / (2. * alpha_sq -  2. * abs((alpha1[nr] / M_PI)) + 1.);
                }
            }

            if(abs(i1 - i2) == 3) {
                if(nr == 0) {
                    alpha2[nr] = alpha_bb_pos * (psi - psi0_bb_pos);
                    alpha3[nr] = alpha_bb_pos;
                    r_min = rmin_pos;
                } else {
                    alpha2[nr] = alpha_bb_neg * (psi - psi0_bb_neg);
                    alpha3[nr] = alpha_bb_neg;
                    r_min = rmin_neg;
                }
            } else {
                alpha2[nr] = alpha_bb_pos * (psi - psi0_bb_pos);
                alpha3[nr] = alpha_bb_pos;
                r_min = rmin_pos;
                if(alpha2[nr] > M_PI || alpha2[nr] < -M_PI) {
                    alpha2[nr] = alpha_bb_neg * (psi - psi0_bb_neg);
                    alpha3[nr] = alpha_bb_neg;
                    r_min = rmin_neg;
                }
            }

            if(alpha2[nr] < M_PI && alpha2[nr] > -M_PI) {
                if(pid_cos) {
                    bb_lambda[nr] = 0.5 * (cos(alpha2[nr]) + 1);
                } else {
                    Real alpha_sq = (alpha2[nr] / M_PI) * (alpha2[nr] / M_PI);
                    bb_lambda[nr] = 1. - alpha_sq 
                                    / (2. * alpha_sq -  2. * abs((alpha2[nr] / M_PI)) + 1.);
                }
            }

            if(bb_lambda[nr] > min_lambda || ss_lambda[nr] > min_lambda) {
                f_var[nr][0] = v4 * v2.norm() /  v4_norm_sq;
                f_var[nr][1] = - v5 * v2.norm() /  v5_norm_sq;
                Real3 df = (f_var[nr][0] * v1.dot(v2) - f_var[nr][1] * v3.dot(v2)) / v2.squaredNorm();
                f_var[nr][2] = - f_var[nr][0] + df;
                f_var[nr][3] = - f_var[nr][1] - df;
            }
        }

        if(norm_too_small)  continue;

        Real ss_lambda_ = ss_lambda[0] * ss_lambda[1];
        Real bb_lambda_ = bb_lambda[0] * bb_lambda[1];

        if(ss_lambda_ < min_lambda && bb_lambda_ < min_lambda)    continue;

        Real force;
        Real lj_energy;

        if(bb_lambda_ > min_lambda) {
            if(dist < r_min * contact_mltp) {
                //TODO update global contact count
            }
            if(sink_pot && dist < r_min * pow(2., 1. / 6.)) {
                lj_energy = - eps_bb;
            }
            else if(pid_barrier) {
                Real rsi = r_min *  pow(2., 1. / 6.) / dist;
                Real r6 = pow(rsi, 6.);
                lj_energy = r6 * (4. * r6 - 18. * rsi + 13.) * eps_bb;
                force += 6. * r6 * (21. * rsi - 8. * r6 - 13.)
                                        / dist * bb_lambda_ * eps_bb;
            } else {
                Real rsi = r_min  / dist;
                Real r6 = pow(rsi, 6.);
                lj_energy = 4. * r6 * (r6 - 1.) * eps_bb;
                force += 24. * r6 * (1. - 2. * r6) / dist * bb_lambda_ * eps_bb;
            }
            //TODO update global potential energy

            for(size_t nr = 0; nr <=1; nr++) {
                size_t other_nr = 1 - nr;
                if(pid_cos) {
                    dvdp[nr] -= 0.5 * alpha3[nr] * sin(alpha2[nr]) * bb_lambda[other_nr] * lj_energy;
                } else {
                    if(alpha2[nr] > 0.) {
                        Real dgdx = 2. * alpha2[nr] * (alpha2[nr] - 1.);
                        Real dgdx2 = (dgdx + 1.) * (dgdx + 1.);
                        dvdp[nr] += alpha3[nr] * dgdx / dgdx2 * bb_lambda[other_nr] * lj_energy;
                    } else {
                        Real dgdx = - 2. * alpha2[nr] * (alpha2[nr] + 1.);
                        Real dgdx2 = (dgdx - 1.) * (dgdx - 1.);
                        dvdp[nr] += alpha3[nr] * dgdx / dgdx2 * bb_lambda[other_nr] * lj_energy;
                    }
                }
            }
        }
        //TODO consider moving analogous pieces of code for bb and ss to external method

        if(eps_mj > min_lambda && ss_lambda_ > min_lambda) {
            //TODO use ssMinimalDistances here for r_min
            if(dist < r_min * contact_mltp) {
                //TODO update global contact count
            }
            //TODO set i_type, j_type, i_charged, j_charged using polarity
            int i_type;
            int j_type;
            bool i_charged;
            bool j_charged;
            if(pid_electrostatics && i_charged && j_charged){
                Real expc = exp(- dist / elektr_screen);
                Real coulpotcoeff = expc * coul;
                if(i_type != j_type)    coulpotcoeff *= -1.;
                if(ele_perm_const) {
                    lj_energy = coulpotcoeff / dist;
                    force += lj_energy * (- dist / elektr_screen - 1.) / dist;
                } else {
                    lj_energy = coulpotcoeff / sq_dist;
                    force += lj_energy * (- dist / elektr_screen - 2.) / dist;
                }
            }
            else if(sink_pot && dist < r_min * pow(2., 1. / 6.)) {
                lj_energy = - eps_mj;
            }
            else if(pid_barrier) {
                Real rsi = r_min *  pow(2., 1. / 6.) / dist;
                Real r6 = pow(rsi, 6.);
                lj_energy = r6 * (4. * r6 - 18. * rsi + 13.) * eps_mj;
                force += 6. * r6 * (21. * rsi - 8. * r6 - 13.) 
                                        / dist * ss_lambda_ * eps_mj;
            } else {
                Real rsi = r_min  / dist;
                Real r6 = pow(rsi, 6.);
                lj_energy = 4. * r6 * (r6 - 1.) * eps_mj;
                force += 24. * r6 * (1. - 2. * r6) / dist * ss_lambda_ * eps_mj;
            }
            //TODO update global potential energy

            for(size_t nr = 0; nr <=1; nr++) {
                size_t other_nr = 1 - nr;
                if(pid_cos) {
                    dvdp[nr] -= 0.5 * alpha_ss* sin(alpha1[nr]) * ss_lambda[other_nr] * lj_energy;
                } else {
                    if(alpha1[nr] > 0.) {
                        Real dgdx = 2. * alpha1[nr] * (alpha1[nr] - 1.);
                        Real dgdx2 = (dgdx + 1.) * (dgdx + 1.);
                        dvdp[nr] += alpha_ss * dgdx / dgdx2 * ss_lambda[other_nr] * lj_energy;
                    } else {
                        Real dgdx = - 2. * alpha1[nr] * (alpha1[nr] + 1.);
                        Real dgdx2 = (dgdx - 1.) * (dgdx - 1.);
                        dvdp[nr] += alpha_ss * dgdx / dgdx2 * ss_lambda[other_nr] * lj_energy;
                    }
                }
            }
        }

        force /= -dist;
        Real3 rep = force * diff_vec;

        forces.col(i) += rep;
        forces.col(j) -= rep;

        for(size_t nr = 0; nr <=1; nr++) {
            int i1 = i;
            int i2 = j;
            if(nr == 1) std::swap(i1,i2);
            forces.col(i1) -= dvdp[nr] * f_var[nr][0];
            forces.col(i1+1) -= dvdp[nr] * f_var[nr][1];
            forces.col(i1-1) -= dvdp[nr] * f_var[nr][2];
            forces.col(i2) -= dvdp[nr] * f_var[nr][3];
        }
    }
    // TODO update energy
}
