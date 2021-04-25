#pragma once
#include <Eigen/Core>

namespace mdk {
    class Topology {
    public:
        Eigen::Vector3d cell, cellInv;

        void setCell(Eigen::Vector3d const& newCell);

        inline Eigen::Vector3d pbc(Eigen::Vector3d v) const {
            for (int d = 0; d < 3; ++d) {
                if (cellInv[d] != 0)
                    v[d] -= cell[d] * (int)(v[d]*cellInv[d]);
            }
            return v;
        }
    };
}
