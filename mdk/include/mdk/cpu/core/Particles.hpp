#pragma once
#include <vector>
#include <mdk/cpu/core/Types.hpp>
#include <mdk/tools/model/CoarseModel.hpp>

namespace mdk {
    class Particles {
    public:
        int n;
        Vectors pos, vel;
        Scalars mass, massInv, charge;
        std::vector<std::pair<int, int>> chains;
        Eigen::Matrix<int8_t, Eigen::Dynamic, 1> types;

        Particles() = default;
        Particles(CoarseModel const& model);
        CoarseModel export() const;
    };
}
