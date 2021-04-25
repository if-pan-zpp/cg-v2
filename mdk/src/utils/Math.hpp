#pragma once
#include <Eigen/Dense>

namespace mdk {
    using VRef = Eigen::Vector3d const&;

    inline double angle(VRef v1, VRef v2, VRef v3) {
        auto u1 = v2 - v1, u2 = v3 - v2;
        return acos(u1.dot(u2)) / (u1.norm() * u2.norm());
    }

    inline double dihedral(VRef v1, VRef v2, VRef v3, VRef v4) {
        auto u1 = v2 - v1, u2 = v3 - v2, u3 = v4 - v3;
        auto d1 = u1.cross(u2), d2 = u3.cross(u2);
        return acos(d1.dot(d2)) / (d1.norm() * d2.norm());
    }
}
