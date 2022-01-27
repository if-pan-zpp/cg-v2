#pragma once
#include <cmath>
#include "utils/Types.hpp"

namespace cg::toolkit {
    inline Real bond(Real3 v1, Real3 v2, Real3 v3) {
        Real3 u1 = v2 - v1, u2 = v3 - v2;
        return acos(-u1.dot(u2) / (u1.norm() * u2.norm()));
    }

    inline Real dihedral(Real3 v1, Real3 v2, Real3 v3, Real3 v4) {
        Real3 u1 = v2 - v1, u2 = v3 - v2, u3 = v4 - v2;
        Real3 d1 = u1.cross(u2), d2 = u2.cross(u3);
        Real phi = acos(d1.dot(d2) / (d1.norm() * d2.norm()));
        if (d1.dot(u3) < 0.) phi = -phi;
        return phi;
    }
}
