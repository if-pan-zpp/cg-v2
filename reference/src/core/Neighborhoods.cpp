#include "core/Geometry.hpp"
#include <stdexcept>

using namespace std;
using namespace CG;

//Neighborhoods Neighborhoods::restrict(Real sub_cutoff, Real sub_pad) const {
//    /* Copy over the same data. */
//    Neighborhoods sub_nb;
//    sub_nb.geometry = geometry;
//    sub_nb.system = geometry->system;
//    sub_nb.cutoff = sub_cutoff;
//    sub_nb.pad = sub_pad;
//    sub_nb.reference_pos = reference_pos;
//    sub_nb.reference_box_shape = reference_box_shape;
//    sub_nb.max_correct_dist = sub_cutoff + 2.0 * sub_pad;
//
//    if (sub_nb.max_correct_dist >= max_correct_dist)
//        throw runtime_error("Geometry - restrict requires strictly lesser top");
//
//    /* As for the pairs, accept only those within the accepted maximal dist. */
//    /* Note: Eigen-based solution would probably be better. */
//    for (auto[p_idx, q_idx]: pairs) {
//        auto p = system->pos.col(p_idx).matrix();
//        auto q = system->pos.col(q_idx).matrix();
//        if (geometry->diff(p, q).norm() <= sub_nb.max_correct_dist) {
//            sub_nb.pairs.emplace_back(p_idx, q_idx);
//        }
//    }
//}

/* Determine the bounding box of a list of positions. */
static pair<Real3, Real3> bbox(Real3List const &ps) {
    auto xmin = ps.row(0).minCoeff();
    auto xmax = ps.row(0).maxCoeff();
    auto ymin = ps.row(1).minCoeff();
    auto ymax = ps.row(1).maxCoeff();
    auto zmin = ps.row(2).minCoeff();
    auto zmax = ps.row(2).maxCoeff();

    auto cmin = Real3(xmin, ymin, zmin);
    auto cmax = Real3(xmax, ymax, zmax);
    return make_pair(cmin, cmax);
}

Neighborhoods::Neighborhoods(Geometry *geometry, Real cutoff, Real pad,
        bool include4) {
    this->cutoff = cutoff;
    this->pad = pad;
    this->geometry = geometry;
    this->include4 = include4;
    this->system = geometry->system;
    reference_pos = system->pos;
    reference_box_shape = geometry->box_shape;
    max_correct_dist = cutoff + 2.0 * pad;

    /* This version of Verlet list is (i) naive, (ii) sort of like in the
     * reference code. We don't update the contacts, though. */
    for (int i = 0; i < system->nresidues; ++i) {
        Real3 pos_i = system->pos.col(i);
        auto cur_chain_id = system->chain_id[i];

        /* Exclude residues on the same chain with abs(j-i) \leq 2. */
        int j;
        for (j = i + 1; j < system->nresidues && j < i + 3 &&
                        system->chain_id[j] == cur_chain_id; ++j);

        for (; j < system->nresidues; ++j) {
            if (j - i == 4 && !include4) continue;

            auto pos_j = system->pos.col(j);
            if (geometry->diff(pos_i, pos_j).norm() <= max_correct_dist)
                pairs.emplace_back(i, j);
        }
    }
}

void Neighborhoods::update_max_correct_dist() {
    /* The formula is:
     * (cutoff+2 pad) - 2 (max particle displacement) - (PBC box change)
     * from what I can see. */
    auto displacements = (reference_pos - system->pos).matrix().colwise();
    Real max_displacement = sqrt(displacements.squaredNorm().maxCoeff());
    Real pbc_box_change = (reference_box_shape - geometry->box_shape).norm();
    max_correct_dist = (cutoff + 2.0*pad)-2.0*max_displacement-pbc_box_change;
    max_correct_dist = max(max_correct_dist, 0.0);
}

void Neighborhoods::update() {
    update_max_correct_dist();
    if (max_correct_dist < cutoff)
        *this = Neighborhoods(this->geometry, cutoff, pad, include4);
}
