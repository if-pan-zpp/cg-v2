#pragma once
#include <eigen3/Eigen/Core>
#include <vector>
#include <string>

namespace cg {
    using Real = double;
    using Real3 = Eigen::Vector3d;
    using Reals3 = Eigen::Matrix3Xd;
    using Reals = Eigen::ArrayXd;

    using Integers = std::vector<int>;
    using Pairs = std::vector<std::pair<int, int>>;
    using Names = std::vector<std::string>;
}
