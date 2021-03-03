#include "core/Geometry.hpp"

using namespace std;
using namespace CG;

Geometry::Geometry(System *system) {
    this->system = system;
    all_neighborhoods = {};

    pbc = {false, false, false};
    box_shape = box_shape_inv = Real3(0, 0, 0);
}

const Neighborhoods *
Geometry::make_neighborhoods(Real cutoff, Real pad, bool include4) {
    auto params = make_tuple(cutoff, pad, include4);

    auto nb_iter = all_neighborhoods.find(params);
    if (nb_iter != all_neighborhoods.end()) {
        return &nb_iter->second;
    } else {
        all_neighborhoods[params] = Neighborhoods(this, cutoff, pad, include4);
        return &all_neighborhoods[params];
    }
}

Real3 Geometry::diff(Real3 p, Real3 q) const {
    Real3 dr = q - p;
    for (int i = 0; i < 3; ++i) {
        if (pbc[i]) {
            /* See the wiki article on periodic boundary condition. */
            dr(i) = abs(dr(i));
            dr(i) -= (int) (dr(i) * box_shape_inv(i) + 0.5) * box_shape(i);
        }
    }
    return dr;
}

void Geometry::set_simulation_box(Real3 box) {
    bool pbc_changed = false;
    for (int i = 0; i < 3; ++i) {
        pbc_changed = pbc_changed || (pbc[i] ^ (box(i) != 0.0));
        box_shape(i) = box(i);
        if (box(i) != 0.0) box_shape_inv(i) = 1.0 / box(i);
    }

    /* If there has been change in the PBC-ness of an axis/axes,
     * recomputation is pretty much necessary. */
    update();
}

void Geometry::update() {
    for (auto& [params, nb]: all_neighborhoods) {
        nb.update();
    }
}